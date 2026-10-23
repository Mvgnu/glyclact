"""Dated continuation of four acquisition paths. Historical runs stay frozen.

Notebook entry point: adaptive_run.run_acquisition_round_20260909.
Only NumPy, pandas and Python's standard library are required.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
from math import erfc, pi, sqrt
from pathlib import Path
import re

import numpy as np
import pandas as pd

import adaptive_run as previous

RUN = "20260909"
LABELS = previous.MEMBER_LABELS
ARMS = ("bayesian_ei", "max_diversity", "bayesian_posterior_mean", "random")
PLATES = {arm: number for number, arm in enumerate(ARMS, 1)}
FEATURES = [f"member_{label}" for label in LABELS]
HERE = Path(__file__).resolve().parent


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_matrix(path):
    """Find exactly 16 contiguous numeric 24-column rows, including temperature exports."""
    rows = list(csv.reader(Path(path).read_text(encoding="utf-8-sig", errors="replace").splitlines()))
    numeric = []
    for index, row in enumerate(rows):
        if len(row) != 24:
            continue
        try:
            values = [float(value.strip()) for value in row]
        except ValueError:
            continue
        if not np.isfinite(values).all():
            raise ValueError(f"Nonfinite reader values in {path}")
        numeric.append((index, values))
    if len(numeric) != 16 or [i for i, _ in numeric] != list(range(numeric[0][0], numeric[0][0] + 16)):
        raise ValueError(f"Expected one contiguous 16x24 plate in {path}")
    wells = [f"{r}{c}" for r in previous.ROWS_384 for c in previous.COLS_384]
    values = np.asarray([v for _, v in numeric]).ravel()
    metadata = "\n".join(row[0] for row in rows[:numeric[0][0]] if row)
    return pd.Series(values, index=wells, name="raw_absorbance"), metadata


def canonical(members):
    labels = members.split("|")
    if len(labels) != len(set(labels)) or not set(labels) <= set(LABELS) or not labels:
        raise ValueError(f"Invalid community members: {members}")
    return previous._members_to_mask(labels)


def load_history(base_dir=HERE, readout_mode="single"):
    """Ordinary mode is available in both rounds; SPR is an explicitly mixed-mode option."""
    base_dir = Path(base_dir)
    if readout_mode not in ("single", "SPR"):
        raise ValueError("readout_mode must be single or SPR")
    old_dir = base_dir / "experiment_setup_output"
    new_dir = old_dir / "20260902_adaptive_12strains"
    old_design_path = old_dir / "community_design_20260828_13strains.csv"
    new_design_path = new_dir / "community_design_all_plates_20260902.csv"
    source_path = new_dir / "source_plate_map_and_dilution_20260902.csv"
    old_readout = base_dir / "readout/Karl_20260830_201548_MAGNUS_ABS-Lactate-assay-rerun.csv"
    source = pd.read_csv(source_path, dtype={"lab_marking": str}).set_index("Source Well").lab_marking.to_dict()
    design2 = pd.read_csv(new_design_path)
    files = [old_design_path, new_design_path, source_path, old_readout]
    datasets = [("20260828", "max_diversity", pd.read_csv(old_design_path), old_readout, "single")]
    # Resolve physical plate IDs from renamed picklists, never an assumed file order.
    for arm, plate in PLATES.items():
        prefix = "maxdiv_repeat" if arm == "max_diversity" else (
            "posterior_mean" if arm == "bayesian_posterior_mean" else arm
        )
        matches = list(new_dir.glob(f"echo_picklist_20260902_{prefix}-plate-*.csv"))
        if len(matches) != 1:
            raise ValueError(f"Cannot uniquely resolve historical picklist for {arm}")
        pick_path = matches[0]
        physical = int(re.search(r"plate-(\d+)", pick_path.name)[1])
        pick = pd.read_csv(pick_path)
        if pick["Destination Plate Name"].nunique() != 1:
            raise ValueError("Historical picklist has multiple destination plates")
        name = pick["Destination Plate Name"].iloc[0]
        sub = design2[design2["Destination Plate Name"].eq(name)].copy()
        for well, g in pick.groupby("Destination Well"):
            condition = sub[sub["Destination Well"].eq(well)]
            if len(condition) != 1 or set(g["Source Well"].map(source)) != set(condition.iloc[0].members.split("|")):
                raise ValueError(f"Picklist/design membership mismatch: {pick_path.name}, {well}")
        if set(pick["Destination Well"]) != set(sub["Destination Well"]):
            raise ValueError("Picklist/design destination mismatch")
        matching_reads = []
        for reader in sorted((base_dir / "readout/04-09").glob("*.csv")):
            _, metadata = read_matrix(reader)
            match = re.search(r"ID1: plate_?(\d+)", metadata)
            mode = "SPR" if "_SPR-" in reader.name else "single"
            if match and int(match[1]) == physical and mode == readout_mode:
                matching_reads.append(reader)
        if len(matching_reads) != 1:
            raise ValueError(f"Expected one {readout_mode} reading for physical plate {physical}")
        reader = matching_reads[0]
        files.extend([pick_path, reader])
        datasets.append(("20260902", arm, sub, reader, readout_mode))

    observations, wells, centers = [], [], []
    for run, arm, design, reader, mode in datasets:
        values, metadata = read_matrix(reader)
        center = float(values.median())
        if design["Destination Well"].duplicated().any() or len(design) != 306:
            raise ValueError(f"Historical layout not 306 unique wells: {run}, {arm}")
        d = design.copy()
        d["community_mask"] = d.members.map(canonical)
        d["raw_absorbance"] = d["Destination Well"].map(values)
        d["run_label"], d["path"], d["readout_mode"] = run, arm, mode
        d["reader_file"], d["reader_sha256"] = reader.name, sha256(reader)
        d["plate_median_384"] = center
        d["response_plate_centered"] = d.raw_absorbance - center
        if d.raw_absorbance.isna().any():
            raise ValueError("Unmatched historical well")
        wells.append(d)
        centers.append({"run_label": run, "path": arm, "reader_file": reader.name,
                        "mode": mode, "plate_median_384": center, "metadata": metadata})
        for mask, g in d.groupby("community_mask"):
            if len(g) != 3 or g.technical_replicate.nunique() != 3:
                raise ValueError("Every historical community must have three replicates")
            observations.append({
                "run_label": run, "path": arm, "community_mask": int(mask),
                "community_id": g.community_id.iloc[0], "members": g.members.iloc[0],
                "community_size": int(g.community_size.iloc[0]),
                "readout_mode": mode, "reader_file": reader.name, "reader_sha256": sha256(reader),
                "destination_wells": "|".join(g["Destination Well"]),
                "raw_replicates": "|".join(f"{v:.6g}" for v in g.raw_absorbance),
                "raw_median": float(g.raw_absorbance.median()), "raw_sd": float(g.raw_absorbance.std()),
                "plate_median_384": center, "response": float(g.response_plate_centered.median()),
            })
    observations = pd.DataFrame(observations)
    if len(observations) != 510:
        raise ValueError(f"Expected 510 historical run/community observations, got {len(observations)}")
    return observations, pd.concat(wells, ignore_index=True), pd.DataFrame(centers), files


def fit_path_models(observations):
    models, training, lineage, summaries = {}, {}, {}, []
    shared = observations[observations.path.eq("max_diversity")]
    for arm in ("bayesian_ei", "bayesian_posterior_mean"):
        rows = pd.concat([shared, observations[observations.path.eq(arm)]], ignore_index=True)
        # A repeated composition contributes a median of its run-level triplicate medians.
        # Keep all original run-level responses in lineage, including disagreement.
        grouped = rows.groupby("community_mask", as_index=False).agg(
            members=("members", "first"), community_size=("community_size", "first"),
            response=("response", "median"), n_observed_runs=("response", "size"),
            min_run_response=("response", "min"), max_run_response=("response", "max"),
            source_runs=("run_label", lambda x: "|".join(x)),
            source_paths=("path", lambda x: "|".join(x)),
        ).sort_values("community_mask").reset_index(drop=True)
        x = previous._mask_matrix(grouped.community_mask.to_numpy(np.uint16))
        for i, feature in enumerate(FEATURES):
            grouped[feature] = x[:, i].astype(int)
        fit, diagnostics = previous.fit_gp(x, grouped.response.to_numpy(float))
        models[arm], training[arm], lineage[arm] = fit, grouped, rows
        summaries.append({
            "model_path": arm, "shared_observations": len(shared),
            "own_acquired_observations": int(rows.path.eq(arm).sum()),
            "total_run_observations": len(rows), "unique_training_compositions": len(grouped),
            "response_definition": "median across runs of triplicate median minus all-384-well plate median",
            "loo_scope": "conditional on hyperparameters selected using all training responses; not nested validation",
            **diagnostics,
        })
        if not set(rows.path) <= {"max_diversity", arm}:
            raise AssertionError("Cross-path outcome leakage")
    return {"models": models, "training": training, "lineage": lineage,
            "model_summary": pd.DataFrame(summaries)}


def stable_ei(mean, sd, incumbent, xi):
    gain = mean - incumbent - xi
    z = gain / sd
    cdf = np.array([0.5 * erfc(-float(v) / sqrt(2)) for v in z])
    return np.maximum(gain * cdf + sd * np.exp(-0.5 * z*z) / sqrt(2*pi), 0)


def available_panel(od_readings, include_296="auto"):
    value = od_readings.get("296")
    finite = value is not None and np.isfinite(float(value))
    if include_296 == "auto":
        if not finite:
            return None
        if float(value) < 0:
            raise ValueError("296 OD cannot be negative; enter 0 for no growth")
        include = float(value) > 0
    elif isinstance(include_296, bool):
        include = include_296
    else:
        raise ValueError("include_296 must be 'auto', True, or False")
    return [label for label in LABELS if label != "296" or include]


def predictions_for_panel(fitted, history, active_labels, xi):
    mask = previous._members_to_mask(active_labels)
    masks = np.arange(1, 1 << len(LABELS), dtype=np.uint16)
    masks = masks[(masks.astype(int) & ~mask) == 0]
    x = previous._mask_matrix(masks)
    candidates = pd.DataFrame({"community_mask": masks.astype(int),
                               "members": ["|".join(previous._mask_to_members(int(m))) for m in masks],
                               "community_size": x.sum(axis=1).astype(int)})
    candidates["previously_picked_any_path"] = candidates.community_mask.isin(set(history.community_mask))
    for arm, short in (("bayesian_ei", "ei"), ("bayesian_posterior_mean", "tm")):
        mean, sd = previous.predict_gp(fitted["models"][arm], x)
        candidates[f"{short}_posterior_mean"] = mean
        candidates[f"{short}_posterior_sd"] = sd
        if short == "ei":
            best = float(fitted["training"][arm].response.max())
            candidates["ei_expected_improvement"] = stable_ei(mean, sd, best, xi)
    candidates["selected_path"] = ""
    candidates["selection_rank"] = pd.Series(pd.NA, index=candidates.index, dtype="Int64")
    return candidates


def maxdiv_continuation(pool, history_masks, count, seed):
    """Balance available richness classes, then continue maximin distance from this path's history."""
    if len(pool) < count:
        raise ValueError("Not enough unseen communities for max diversity")
    p = pool.copy().reset_index(drop=True)
    x = previous._mask_matrix(p.community_mask.to_numpy(np.uint16)).astype(np.int16)
    hist = previous._mask_matrix(np.array(sorted(set(history_masks)), dtype=np.uint16)).astype(np.int16)
    nearest = (x.sum(1)[:, None] + hist.sum(1)[None, :] - 2*x@hist.T).min(axis=1)
    member_counts = hist.sum(axis=0)
    sizes = p.community_size.to_numpy()
    capacities = p.community_size.value_counts().to_dict()
    quota = {k: 0 for k in sorted(capacities)}
    for _ in range(count):
        k = min((k for k in quota if quota[k] < capacities[k]), key=lambda k: (quota[k], k))
        quota[k] += 1
    rng = np.random.default_rng(seed)
    selected = []
    used = np.zeros(len(p), dtype=bool)
    for k, take in quota.items():
        for _ in range(take):
            eligible = np.flatnonzero((sizes == k) & ~used)
            eligible = eligible[nearest[eligible] == nearest[eligible].max()]
            reuse = x[eligible] @ member_counts
            eligible = eligible[reuse == reuse.min()]
            choice = int(rng.choice(eligible))
            selected.append(choice)
            used[choice] = True
            member_counts += x[choice]
            nearest = np.minimum(nearest, np.abs(x-x[choice]).sum(axis=1))
    return p.iloc[selected].copy().reset_index(drop=True)


def select_paths(candidates, history, seed=20260909, count=102):
    """Continue previous dedup priority: maxdiv -> EI -> TM -> random; retain pre-dedup proposals."""
    pool = candidates[~candidates.previously_picked_any_path].copy()
    maxdiv = maxdiv_continuation(pool, history.loc[history.path.eq("max_diversity"), "community_mask"], count, seed)
    rankings = {
        "bayesian_ei": pool.sort_values(["ei_expected_improvement", "ei_posterior_mean", "community_mask"], ascending=[False, False, True]),
        "bayesian_posterior_mean": pool.sort_values(["tm_posterior_mean", "tm_posterior_sd", "community_mask"], ascending=[False, True, True]),
        "random": pool.iloc[np.random.default_rng(seed + 1).permutation(len(pool))],
    }
    chosen = {"max_diversity": maxdiv}
    used = set(maxdiv.community_mask)
    proposal_rows = []
    for arm in ("bayesian_ei", "bayesian_posterior_mean", "random"):
        ranking = rankings[arm]
        picked = ranking[~ranking.community_mask.isin(used)].head(count).copy()
        for rank, row in enumerate(ranking.head(count).itertuples(), 1):
            proposal_rows.append({"path": arm, "pre_dedup_rank": rank,
                                  "community_mask": row.community_mask,
                                  "blocked_by_earlier_current_path": row.community_mask in used,
                                  "retained_in_this_path": row.community_mask in set(picked.community_mask)})
        if len(picked) != count:
            raise ValueError(f"Not enough candidates for {arm}")
        chosen[arm] = picked.reset_index(drop=True)
        used.update(picked.community_mask)
    result = candidates.copy()
    for arm, batch in chosen.items():
        batch["selected_path"] = arm
        batch["selection_rank"] = range(1, len(batch) + 1)
        ranks = batch.set_index("community_mask").selection_rank
        match = result.community_mask.isin(ranks.index)
        result.loc[match, "selected_path"] = arm
        result.loc[match, "selection_rank"] = result.loc[match, "community_mask"].map(ranks).astype("Int64")
    assert len(used) == 4*count
    assert used.isdisjoint(set(history.community_mask))
    return chosen, result, pd.DataFrame(proposal_rows)


def layout(chosen, active_labels, seed):
    parts = []
    for arm in ARMS:
        number = PLATES[arm]
        p = previous._layout_new_batch(chosen[arm], arm,
            f"community_{RUN}_{arm}_plate_{number}", seed + 100 + number, 3)
        # Preserve the historical function's field order for Echo's transfer helper.
        p["community_id"] = p.community_id.map(lambda cid: f"R3_{cid}")
        p["repeat_status"] = "new_composition"
        parts.append(p)
    design = pd.concat(parts, ignore_index=True)
    for label in LABELS:
        design[f"member_{label}"] = design.members.str.split("|").map(lambda xs: int(label in xs))
    assert len(design) == 1224 and design.community_mask.nunique() == 408
    assert not design.duplicated(["Destination Plate Name", "Destination Well"]).any()
    assert set(design["Destination Well"]).isdisjoint(previous.EDGE_WELLS)
    assert design.groupby("community_mask").technical_replicate.nunique().eq(3).all()
    assert design.members.str.split("|").map(lambda xs: set(xs) <= set(active_labels)).all()
    return design


def serializable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def make_fixed_dose_picklists(design, source, source_plate_name, source_load_ul):
    """Keep September 2's executed dose; no 15 uL rule or automatic biomass switch."""
    invalid = source.loc[~np.isfinite(source.working_source_OD600)
                         | source.working_source_OD600.le(0), "lab_marking"].tolist()
    if invalid:
        return False, "OD input required for: " + ", ".join(invalid), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    transfers = previous._transfer_rows(
        design, source, source_plate_name, 0.0, 0.100, "fixed_per_member_biomass", 25)
    audit = transfers.groupby(
        ["Source Strategy", "Destination Plate Name", "Destination Well", "community_id",
         "technical_replicate", "community_size", "dosing_mode"], as_index=False,
    ).agg(total_inoculum_nL=("executed_transfer_nL", "sum"),
          realized_total_OD600_uL=("realized_member_OD600_uL", "sum"))
    audit["target_total_OD600_uL"] = 0.100 * audit.community_size
    audit["biomass_error_pct"] = 100 * (audit.realized_total_OD600_uL / audit.target_total_OD600_uL - 1)
    demand = transfers.groupby(
        ["Destination Plate Name", "lab_marking", "Source Well"], as_index=False,
    ).agg(destination_wells=("Destination Well", "nunique"),
          transfer_demand_nL=("executed_transfer_nL", "sum"))
    demand["transfer_demand_uL"] = demand.transfer_demand_nL / 1000
    demand["source_load_uL"] = source_load_ul
    demand["expected_remaining_uL"] = source_load_ul - demand.transfer_demand_uL
    # Bookkeeping, not an instrument dead-volume specification. Never silently change dose.
    demand["within_loaded_volume"] = demand.expected_remaining_uL.ge(0)
    if transfers.executed_transfer_nL.le(0).any():
        return False, "At least one calculated transfer rounded to zero.", transfers, audit, demand
    if not demand.within_loaded_volume.all():
        return False, "Source demand exceeds the loaded liquid; dose unchanged. Resolve source loading before export.", transfers, audit, demand
    message = (
        "Picklists ready: fixed 0.100 OD600·µL per member, rounded to 25 nL; "
        "no 15 µL cap/reserve or automatic dosing fallback. "
        "Restore the source row to its configured fill before each destination plate. "
        "Remaining volumes are bookkeeping estimates, not confirmation of the source plate's ejectable range."
    )
    return True, message, transfers, audit, demand


def run_round(*, od_readings, include_296="auto", readout_mode="single", seed=20260909,
              source_plate_name="396", source_load_ul=50.0, target_working_od=1.2,
              medium_ul=50.0,
              xi=0.002, base_dir=HERE, export=True, fitted_history=None):
    """Train first; pending 296 produces only labelled scenario designs, never an Echo file."""
    base_dir = Path(base_dir)
    if set(od_readings) - set(LABELS):
        raise ValueError("OD dictionary contains an unknown species")
    if not (0 < source_load_ul <= 50):
        raise ValueError("Require 0 < source fill <= 50 µL")
    if not np.isfinite([target_working_od, medium_ul]).all() or target_working_od <= 0 or medium_ul <= 0:
        raise ValueError("Volumes and target working OD must be positive")
    active = available_panel(od_readings, include_296)
    if fitted_history is None:
        observations, wells, centers, files = load_history(base_dir, readout_mode)
        fitted = fit_path_models(observations)
    else:
        observations, wells, centers, files, fitted = fitted_history
        modes = set(observations.loc[observations.run_label.eq("20260902"), "readout_mode"])
        if modes != {readout_mode}:
            raise ValueError("Cached history readout mode mismatch")
    panels = [active] if active is not None else [[s for s in LABELS if s != "296"], LABELS.copy()]
    scenarios = {}
    for labels in panels:
        candidates = predictions_for_panel(fitted, observations, labels, xi)
        chosen, candidates, overlap = select_paths(candidates, observations, seed)
        design = layout(chosen, labels, seed)
        source = previous.build_source_plan(od_readings, labels, source_load_ul, 10.0, target_working_od)
        ready, message, transfers, audit, demand = make_fixed_dose_picklists(
            design, source, source_plate_name, source_load_ul)
        if active is None:
            ready, message = False, "296 OD pending: scenario preview only; enter 0 to exclude or a positive reading to include."
        if len(audit):
            audit["media_prefill_uL"] = medium_ul
            audit["final_culture_volume_uL"] = medium_ul + audit.total_inoculum_nL/1000
        if len(transfers):
            assert transfers["Transfer Volume"].mod(25).eq(0).all()
            assert set(transfers["lab_marking"]) <= set(labels)
        summary = design.groupby(["Source Strategy", "Destination Plate Name"], as_index=False).agg(
            communities=("community_mask", "nunique"), culture_wells=("Destination Well", "nunique"),
            min_richness=("community_size", "min"), max_richness=("community_size", "max"))
        scenarios[f"{len(labels)}strains"] = {"active_labels": labels, "candidates": candidates,
            "design": design, "source_plan": source, "overlap_audit": overlap,
            "summary": summary, "picklist_ready": bool(ready), "status": message,
            "transfers": transfers, "well_audit": audit, "source_demand": demand}

    result = {"model_summary": fitted["model_summary"], "scenarios": scenarios,
              "fitted_history": (observations, wells, centers, files, fitted), "output_dir": None}
    if not export:
        return result
    # Each execution is a new immutable snapshot, including pending-OD scenarios.
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    output = base_dir / "experiment_setup_output" / f"{RUN}_paths" / stamp
    output.mkdir(parents=True, exist_ok=False)
    result["output_dir"] = output
    observations.to_csv(output / "history_run_observations.csv", index=False)
    wells.to_csv(output / "history_wells.csv", index=False)
    centers.to_csv(output / "readout_offsets_and_metadata.csv", index=False)
    fitted["model_summary"].to_csv(output / "model_summary.csv", index=False)
    for arm in fitted["models"]:
        fitted["training"][arm].to_csv(output / f"training_{arm}.csv", index=False)
        fitted["lineage"][arm].to_csv(output / f"training_lineage_{arm}.csv", index=False)
        model = fitted["models"][arm]
        np.savez_compressed(output / f"model_{arm}.npz", x_train=model.x_train,
            y_standardized=model.y_train_standardized, chol=model.chol, alpha=model.alpha,
            feature_order=np.asarray(LABELS), y_mean=model.y_mean, y_scale=model.y_scale,
            length_scale=model.length_scale, signal_variance=model.signal_variance,
            noise_variance=model.noise_variance)
    for scenario, data in scenarios.items():
        subdir = output / scenario
        subdir.mkdir()
        for key in ("candidates", "design", "source_plan", "overlap_audit", "summary", "well_audit", "source_demand"):
            if len(data[key]):
                data[key].to_csv(subdir / f"{key}.csv", index=False)
        if len(data["transfers"]):
            data["transfers"].to_csv(subdir / "transfer_audit.csv", index=False)
        for arm in ARMS:
            data["design"][data["design"]["Source Strategy"].eq(arm)].to_csv(
                subdir / f"community_design_{RUN}_{arm}-plate-{PLATES[arm]}.csv", index=False)
        if data["picklist_ready"]:
            columns = ["Source Plate Name", "Source Well", "Destination Plate Name", "Destination Well", "Transfer Volume"]
            for arm in ARMS:
                data["transfers"][data["transfers"]["Source Strategy"].eq(arm)][columns].to_csv(
                    subdir / f"echo_picklist_{RUN}_{arm}-plate-{PLATES[arm]}.csv", index=False)
    # Snapshot executable code plus complete input hashes, not only the resulting CSVs.
    for code in (Path(__file__), base_dir / "adaptive_run.py"):
        (output / code.name).write_bytes(code.read_bytes())
    manifest = {
        "run_label": RUN, "snapshot_utc": stamp,
        "od_readings_at_1_to_10": {s: (float(od_readings[s]) if od_readings.get(s) is not None and np.isfinite(float(od_readings[s])) else None) for s in LABELS},
        "include_296": include_296, "readout_mode": readout_mode,
        "centering": "all_384_well_plate_median", "seed": seed, "ei_xi": xi,
        "source_plate_name": source_plate_name, "source_load_ul": source_load_ul,
        "source_dead_volume_ul": None, "target_working_od": target_working_od,
        "max_inoculum_ul": None, "medium_prefill_ul": medium_ul,
        "dosing_mode": "fixed_per_member_biomass", "member_biomass_OD600_uL": 0.100,
        "automatic_dosing_fallback": False,
        "technical_replicates": 3, "quota_per_path": 102,
        "shared_data": "both max-diversity rounds; median across repeated compositions",
        "dedup_priority": ["max_diversity", "bayesian_ei", "bayesian_posterior_mean", "random"],
        "input_files": [{"path": str(p.relative_to(base_dir)), "sha256": sha256(p)} for p in files],
        "code_hashes": {p.name: sha256(p) for p in (Path(__file__), base_dir / "adaptive_run.py")},
        "numpy_version": np.__version__, "pandas_version": pd.__version__,
        "scenarios": {name: {"ready": d["picklist_ready"], "status": d["status"], "active_species": d["active_labels"]} for name, d in scenarios.items()},
    }
    manifest["artifact_sha256"] = {str(p.relative_to(output)): sha256(p) for p in output.rglob("*") if p.is_file()}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, default=serializable, allow_nan=False) + "\n")
    return result


if __name__ == "__main__":
    result = run_round(od_readings={label: None for label in LABELS})
    print(result["model_summary"].to_string(index=False))
    for name, data in result["scenarios"].items():
        print(name, data["status"])
        print(data["summary"].to_string(index=False))
    print(result["output_dir"])
