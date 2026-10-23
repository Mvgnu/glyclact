"""Adaptive four-plate Echo design used by experiment_setup.ipynb.

The module intentionally uses only NumPy and pandas so the runbook does not depend on
PyMC/scikit-learn being installed on the lab computer.  The Bayesian surrogate is a
Gaussian process over the 13-member binary composition vector with an RBF/Hamming
kernel.  Hyperparameters are selected by marginal likelihood.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import erf, exp, log, pi, sqrt
from pathlib import Path
import itertools
import string

import numpy as np
import pandas as pd


MEMBER_LABELS = [
    "256", "75", "269", "294", "275", "41", "158",
    "145", "200", "196", "291", "296", "264",
]

ROWS_384 = list(string.ascii_uppercase[:16])
COLS_384 = list(range(1, 25))
INTERIOR_WELLS = [f"{row}{col}" for row in ROWS_384[1:-1] for col in COLS_384[1:-1]]
EDGE_WELLS = {
    f"{row}{col}"
    for row in ROWS_384
    for col in COLS_384
    if row in ("A", "P") or col in (1, 24)
}
# The six assay standards occupy the left column, A1:F1.  All other edge wells
# are same-medium zero controls.
STANDARD_WELLS = {f"{row}1" for row in ROWS_384[:6]}


@dataclass
class GPFit:
    length_scale: float
    signal_variance: float
    noise_variance: float
    log_marginal_likelihood: float
    y_mean: float
    y_scale: float
    x_train: np.ndarray
    y_train_standardized: np.ndarray
    chol: np.ndarray
    alpha: np.ndarray


def _members_to_mask(members: str | tuple[str, ...] | list[str]) -> int:
    labels = members.split("|") if isinstance(members, str) else list(members)
    label_set = set(labels)
    return sum((1 << i) for i, label in enumerate(MEMBER_LABELS) if label in label_set)


def _mask_to_members(mask: int) -> tuple[str, ...]:
    return tuple(label for i, label in enumerate(MEMBER_LABELS) if mask & (1 << i))


def _mask_matrix(masks: np.ndarray) -> np.ndarray:
    shifts = np.arange(len(MEMBER_LABELS), dtype=np.uint16)
    return ((masks[:, None].astype(np.uint16) >> shifts) & 1).astype(float)


def _load_reader_plate(path: Path) -> dict[str, float]:
    raw = pd.read_csv(path, header=None, skiprows=9, nrows=16, engine="python")
    if raw.shape[1] < 24:
        raise ValueError(f"Reader file has {raw.shape[1]} columns; expected 24: {path}")
    values = raw.iloc[:, :24].apply(pd.to_numeric, errors="coerce")
    if values.isna().any().any():
        raise ValueError(f"Reader matrix contains non-numeric/missing values: {path}")
    return {
        f"{row}{col}": float(values.iloc[row_i, col - 1])
        for row_i, row in enumerate(ROWS_384)
        for col in COLS_384
    }


def load_training_data(design_path: Path, readout_path: Path) -> tuple[pd.DataFrame, float]:
    """Join the previous design to the rerun and robustly aggregate triplicates."""
    design = pd.read_csv(design_path, dtype={"lab_marking": str})
    plate = _load_reader_plate(readout_path)
    design["raw_absorbance"] = design["Destination Well"].map(plate)
    if design["raw_absorbance"].isna().any():
        missing = design.loc[design["raw_absorbance"].isna(), "Destination Well"].tolist()
        raise ValueError(f"Readout missing designed wells: {missing[:10]}")

    blank_wells = sorted(EDGE_WELLS - STANDARD_WELLS)
    blank_median = float(np.median([plate[well] for well in blank_wells]))

    grouped = design.groupby(
        ["community_id", "community_size", "members"], as_index=False
    )["raw_absorbance"].agg(
        response_raw_median="median",
        replicate_mean="mean",
        replicate_min="min",
        replicate_max="max",
        replicate_sd="std",
        n_replicates="count",
    )
    grouped["response_blank_corrected"] = grouped["response_raw_median"] - blank_median
    grouped["replicate_range"] = grouped["replicate_max"] - grouped["replicate_min"]
    grouped["community_mask"] = grouped["members"].map(_members_to_mask)
    for label in MEMBER_LABELS:
        grouped[f"member_{label}"] = grouped["members"].str.split("|").map(
            lambda values, label=label: int(label in values)
        )

    if not grouped["community_mask"].is_unique:
        raise ValueError("Previous design contains duplicate community compositions.")
    if not grouped["n_replicates"].eq(3).all():
        raise ValueError("Expected exactly three technical replicates per community.")
    return grouped.sort_values("community_mask").reset_index(drop=True), blank_median


def _rbf_hamming(x_a: np.ndarray, x_b: np.ndarray, length_scale: float) -> np.ndarray:
    # Squared Euclidean distance between binary vectors is Hamming distance.
    distance = (
        np.sum(x_a * x_a, axis=1)[:, None]
        + np.sum(x_b * x_b, axis=1)[None, :]
        - 2.0 * x_a @ x_b.T
    )
    return np.exp(-np.maximum(distance, 0.0) / (2.0 * length_scale**2))


def fit_gp(x_train: np.ndarray, y_train: np.ndarray) -> tuple[GPFit, dict[str, float]]:
    """Fit a small exact GP by marginal-likelihood grid search."""
    y_mean = float(np.mean(y_train))
    y_scale = float(np.std(y_train, ddof=1))
    if not np.isfinite(y_scale) or y_scale <= 0:
        y_scale = 1.0
    y = (y_train - y_mean) / y_scale

    length_scales = [0.35, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00]
    signal_variances = [0.50, 1.00, 2.00]
    noise_variances = [0.01, 0.03, 0.10, 0.30, 0.60]
    best = None
    n = len(y)
    eye = np.eye(n)

    for length_scale, signal_variance, noise_variance in itertools.product(
        length_scales, signal_variances, noise_variances
    ):
        kernel = signal_variance * _rbf_hamming(x_train, x_train, length_scale)
        covariance = kernel + (noise_variance + 1e-8) * eye
        try:
            chol = np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError:
            continue
        alpha = np.linalg.solve(chol.T, np.linalg.solve(chol, y))
        log_likelihood = float(
            -0.5 * y @ alpha
            - np.sum(np.log(np.diag(chol)))
            - 0.5 * n * log(2.0 * pi)
        )
        if best is None or log_likelihood > best[0]:
            best = (
                log_likelihood, length_scale, signal_variance,
                noise_variance, chol, alpha,
            )

    if best is None:
        raise RuntimeError("No numerically valid GP hyperparameter combination was found.")

    log_likelihood, length_scale, signal_variance, noise_variance, chol, alpha = best
    fit = GPFit(
        length_scale=length_scale,
        signal_variance=signal_variance,
        noise_variance=noise_variance,
        log_marginal_likelihood=log_likelihood,
        y_mean=y_mean,
        y_scale=y_scale,
        x_train=x_train,
        y_train_standardized=y,
        chol=chol,
        alpha=alpha,
    )

    # Exact leave-one-out predictions for a Gaussian model.
    covariance_inverse = np.linalg.solve(chol.T, np.linalg.solve(chol, eye))
    loo_standardized = y - alpha / np.diag(covariance_inverse)
    loo = y_mean + y_scale * loo_standardized
    residual = y_train - loo
    pearson = float(pd.Series(y_train).corr(pd.Series(loo), method="pearson"))
    # Spearman is Pearson correlation of ranks; spelling it out avoids a SciPy dependency.
    spearman = float(
        pd.Series(y_train).rank(method="average").corr(
            pd.Series(loo).rank(method="average"), method="pearson"
        )
    )
    diagnostics = {
        "n_training_communities": float(n),
        "length_scale": length_scale,
        "signal_variance_standardized": signal_variance,
        "noise_variance_standardized": noise_variance,
        "log_marginal_likelihood": log_likelihood,
        "loo_rmse_absorbance": float(sqrt(np.mean(residual**2))),
        "loo_mae_absorbance": float(np.mean(np.abs(residual))),
        "loo_pearson": pearson,
        "loo_spearman": spearman,
    }
    return fit, diagnostics


def predict_gp(fit: GPFit, x_candidate: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    cross_kernel = fit.signal_variance * _rbf_hamming(
        x_candidate, fit.x_train, fit.length_scale
    )
    mean_standardized = cross_kernel @ fit.alpha
    solved = np.linalg.solve(fit.chol, cross_kernel.T)
    variance_standardized = np.maximum(
        fit.signal_variance - np.sum(solved * solved, axis=0), 1e-12
    )
    mean = fit.y_mean + fit.y_scale * mean_standardized
    sd = fit.y_scale * np.sqrt(variance_standardized)
    return mean, sd


def _normal_cdf(values: np.ndarray) -> np.ndarray:
    return np.fromiter(
        (0.5 * (1.0 + erf(float(value) / sqrt(2.0))) for value in values),
        dtype=float,
        count=len(values),
    )


def expected_improvement(
    mean: np.ndarray,
    sd: np.ndarray,
    best_observed: float,
    xi: float = 0.002,
) -> np.ndarray:
    improvement = mean - best_observed - xi
    z = np.divide(improvement, sd, out=np.zeros_like(improvement), where=sd > 0)
    cdf = _normal_cdf(z)
    pdf = np.exp(-0.5 * z**2) / sqrt(2.0 * pi)
    result = improvement * cdf + sd * pdf
    result[sd <= 0] = np.maximum(improvement[sd <= 0], 0.0)
    return np.maximum(result, 0.0)


def build_candidate_predictions(
    training: pd.DataFrame,
    fit: GPFit,
    ei_xi: float,
    active_labels: list[str],
) -> pd.DataFrame:
    active_mask = _members_to_mask(active_labels)
    masks = np.arange(1, 1 << len(MEMBER_LABELS), dtype=np.uint16)
    masks = masks[(masks.astype(int) & ~active_mask) == 0]
    matrix = _mask_matrix(masks)
    mean, sd = predict_gp(fit, matrix)
    ei = expected_improvement(
        mean, sd, float(training["response_blank_corrected"].max()), xi=ei_xi
    )
    predictions = pd.DataFrame({
        "community_mask": masks.astype(int),
        "community_size": matrix.sum(axis=1).astype(int),
        "members": ["|".join(_mask_to_members(int(mask))) for mask in masks],
        "posterior_mean_blank_corrected": mean,
        "posterior_sd": sd,
        "expected_improvement": ei,
    })
    tested_masks = set(training["community_mask"].astype(int))
    predictions["tested_in_previous_run"] = predictions["community_mask"].isin(tested_masks)
    predictions["posterior_mean_rank"] = predictions[
        "posterior_mean_blank_corrected"
    ].rank(method="first", ascending=False).astype(int)
    predictions["expected_improvement_rank"] = predictions[
        "expected_improvement"
    ].rank(method="first", ascending=False).astype(int)
    predictions["selected_strategy"] = ""
    predictions["selection_rank"] = pd.Series(pd.NA, index=predictions.index, dtype="Int64")
    return predictions


def select_strategies(
    predictions: pd.DataFrame,
    previous_design: pd.DataFrame,
    n_per_strategy: int,
    random_seed: int,
    active_labels: list[str],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Select four disjoint batches; the repeat is the only historical reuse.

    Deduplication priority is repeat -> EI -> posterior mean -> random.  Model-based
    strategies therefore retain their highest-ranked unique proposals and random cannot
    accidentally steal them first.
    """
    previous_unique = (
        previous_design[["community_id", "community_size", "members"]]
        .drop_duplicates("community_id")
        .copy()
    )
    previous_unique["community_mask"] = previous_unique["members"].map(_members_to_mask)
    if len(previous_unique) != n_per_strategy:
        raise ValueError(
            f"Expected {n_per_strategy} prior communities for exact repeat; found {len(previous_unique)}."
        )

    active_mask = _members_to_mask(active_labels)
    previous_unique["is_feasible_exact_repeat"] = previous_unique["community_mask"].map(
        lambda mask: (int(mask) & ~active_mask) == 0
    )
    selected_masks = set(
        previous_unique.loc[previous_unique["is_feasible_exact_repeat"], "community_mask"].astype(int)
    )

    # Keep every feasible historical community as an exact repeat.  Each slot that
    # previously contained an unavailable member gets a same-richness, unseen,
    # maximin-Hamming backfill over the active panel.
    repeat_rows = []
    member_counts = _mask_matrix(np.array(sorted(selected_masks), dtype=np.uint16)).sum(axis=0)
    rng = np.random.default_rng(random_seed)
    for slot in previous_unique.itertuples(index=False):
        if slot.is_feasible_exact_repeat:
            prediction = predictions.loc[predictions["community_mask"].eq(int(slot.community_mask))].iloc[0].to_dict()
            prediction.update({
                "community_id": slot.community_id,
                "slot_historical_community_id": slot.community_id,
                "repeat_status": "exact_repeat",
            })
            repeat_rows.append(prediction)
            continue

        target_size = min(int(slot.community_size), len(active_labels))
        pool = predictions.loc[
            predictions["community_size"].eq(target_size)
            & ~predictions["tested_in_previous_run"]
            & ~predictions["community_mask"].isin(selected_masks)
        ].copy()
        if pool.empty:
            pool = predictions.loc[
                ~predictions["tested_in_previous_run"]
                & ~predictions["community_mask"].isin(selected_masks)
            ].copy()
        pool_matrix = _mask_matrix(pool["community_mask"].to_numpy(np.uint16))
        if selected_masks:
            selected_matrix = _mask_matrix(np.array(sorted(selected_masks), dtype=np.uint16))
            distances = (
                pool_matrix[:, None, :] != selected_matrix[None, :, :]
            ).sum(axis=2)
            nearest_distance = distances.min(axis=1)
        else:
            nearest_distance = np.full(len(pool), len(active_labels))
        eligible_positions = np.flatnonzero(nearest_distance == nearest_distance.max())
        reuse_scores = pool_matrix[eligible_positions] @ member_counts
        eligible_positions = eligible_positions[reuse_scores == reuse_scores.min()]
        chosen_position = int(rng.choice(eligible_positions))
        chosen = pool.iloc[chosen_position].to_dict()
        chosen_mask = int(chosen["community_mask"])
        selected_masks.add(chosen_mask)
        member_counts += pool_matrix[chosen_position]
        chosen.update({
            "community_id": f"MAXDIV_BACKFILL_{len(repeat_rows) + 1:03d}",
            "slot_historical_community_id": slot.community_id,
            "repeat_status": "active_panel_maxdiv_backfill",
        })
        repeat_rows.append(chosen)

    repeat_batch = pd.DataFrame(repeat_rows)
    if len(repeat_batch) != n_per_strategy or not repeat_batch["community_mask"].is_unique:
        raise AssertionError("Repeat/backfill plate is not a full unique batch.")

    eligible = predictions.loc[
        ~predictions["tested_in_previous_run"]
        & ~predictions["community_mask"].isin(selected_masks)
    ].copy()

    ei = eligible.sort_values(
        ["expected_improvement", "posterior_mean_blank_corrected", "community_mask"],
        ascending=[False, False, True],
    ).head(n_per_strategy).copy()
    selected_masks.update(ei["community_mask"].astype(int))

    posterior_mean = eligible.loc[
        ~eligible["community_mask"].isin(selected_masks)
    ].sort_values(
        ["posterior_mean_blank_corrected", "posterior_sd", "community_mask"],
        ascending=[False, True, True],
    ).head(n_per_strategy).copy()
    selected_masks.update(posterior_mean["community_mask"].astype(int))

    random_pool = eligible.loc[~eligible["community_mask"].isin(selected_masks)].copy()
    chosen_index = rng.choice(random_pool.index.to_numpy(), size=n_per_strategy, replace=False)
    random_batch = random_pool.loc[chosen_index].copy()

    batches = {
        "max_diversity_repeat": repeat_batch.reset_index(drop=True),
        "bayesian_ei": ei.reset_index(drop=True),
        "bayesian_posterior_mean": posterior_mean.reset_index(drop=True),
        "random": random_batch.reset_index(drop=True),
    }

    for strategy, batch in batches.items():
        if strategy == "max_diversity_repeat":
            masks = batch["community_mask"].astype(int).tolist()
            ranks = list(range(1, len(batch) + 1))
        else:
            sort_column = {
                "bayesian_ei": "expected_improvement",
                "bayesian_posterior_mean": "posterior_mean_blank_corrected",
                "random": "community_mask",
            }[strategy]
            ascending = strategy == "random"
            batch.sort_values(sort_column, ascending=ascending, inplace=True)
            batch.reset_index(drop=True, inplace=True)
            masks = batch["community_mask"].astype(int).tolist()
            ranks = list(range(1, len(batch) + 1))
        mask_to_rank = dict(zip(masks, ranks))
        predictions.loc[predictions["community_mask"].isin(masks), "selected_strategy"] = strategy
        predictions.loc[predictions["community_mask"].isin(masks), "selection_rank"] = (
            predictions.loc[predictions["community_mask"].isin(masks), "community_mask"]
            .map(mask_to_rank)
            .astype("Int64")
        )

    all_current_masks = [
        int(mask)
        for batch in batches.values()
        for mask in batch["community_mask"].tolist()
    ]
    if len(all_current_masks) != len(set(all_current_masks)):
        raise AssertionError("A community was selected by more than one current strategy.")
    return batches, predictions


def _layout_new_batch(
    batch: pd.DataFrame,
    strategy: str,
    plate_name: str,
    random_seed: int,
    technical_replicates: int,
) -> pd.DataFrame:
    n_conditions = len(batch)
    if n_conditions * technical_replicates > len(INTERIOR_WELLS):
        raise ValueError("Batch does not fit in the interior destination wells.")
    rng = np.random.default_rng(random_seed)
    condition_order = rng.permutation(n_conditions)
    offsets = [i * (n_conditions + 1) for i in range(technical_replicates)]
    rows = []
    for slot, condition_i in enumerate(condition_order):
        condition = batch.iloc[condition_i]
        replicate_labels = rng.permutation(np.arange(1, technical_replicates + 1))
        for block_i, offset in enumerate(offsets):
            rows.append({
                "Source Strategy": strategy,
                "repeat_status": "not_applicable",
                "Destination Plate Name": plate_name,
                "Destination Well": INTERIOR_WELLS[slot + offset],
                "technical_replicate": int(replicate_labels[block_i]),
                "community_id": f"{strategy.upper()}_{condition_i + 1:03d}",
                "community_mask": int(condition["community_mask"]),
                "community_size": int(condition["community_size"]),
                "members": condition["members"],
            })
    return pd.DataFrame(rows)


def build_plate_designs(
    batches: dict[str, pd.DataFrame],
    previous_design: pd.DataFrame,
    run_label: str,
    random_seed: int,
    technical_replicates: int = 3,
) -> pd.DataFrame:
    designs = []
    repeat_rows = []
    for condition in batches["max_diversity_repeat"].itertuples(index=False):
        slot_rows = previous_design.loc[
            previous_design["community_id"].eq(condition.slot_historical_community_id)
        ]
        if len(slot_rows) != technical_replicates:
            raise AssertionError("Historical repeat slot does not contain three replicates.")
        for slot in slot_rows.itertuples(index=False):
            repeat_rows.append({
                "Source Strategy": "max_diversity_repeat",
                "repeat_status": condition.repeat_status,
                "Destination Plate Name": f"community_{run_label}_maxdiv_repeat",
                "Destination Well": getattr(slot, "_1"),
                "technical_replicate": slot.technical_replicate,
                "community_id": condition.community_id,
                "community_mask": int(condition.community_mask),
                "community_size": int(condition.community_size),
                "members": condition.members,
            })
    designs.append(pd.DataFrame(repeat_rows))

    plate_suffixes = {
        "bayesian_ei": "bayesian_ei",
        "bayesian_posterior_mean": "posterior_mean",
        "random": "random",
    }
    for offset, strategy in enumerate(("bayesian_ei", "bayesian_posterior_mean", "random"), 1):
        designs.append(_layout_new_batch(
            batches[strategy],
            strategy,
            f"community_{run_label}_{plate_suffixes[strategy]}",
            random_seed + offset,
            technical_replicates,
        ))

    design = pd.concat(designs, ignore_index=True, sort=False)
    base_columns = [
        "Source Strategy", "repeat_status", "Destination Plate Name", "Destination Well",
        "technical_replicate", "community_id", "community_mask",
        "community_size", "members",
    ]
    design = design[base_columns]
    for label in MEMBER_LABELS:
        design[f"member_{label}"] = design["members"].str.split("|").map(
            lambda values, label=label: int(label in values)
        )

    if set(design["Destination Well"]) & EDGE_WELLS:
        raise AssertionError("A design uses edge wells.")
    well_counts = design.groupby("Destination Plate Name")["Destination Well"].agg(
        rows="size", unique_wells="nunique"
    )
    if not well_counts["rows"].eq(well_counts["unique_wells"]).all():
        raise AssertionError("Duplicate destination wells within a plate.")
    if not design.groupby(["Destination Plate Name", "community_mask"]).size().eq(3).all():
        raise AssertionError("Every community must have three technical replicates.")
    per_strategy_masks = design.groupby("Source Strategy")["community_mask"].unique()
    flattened = [int(mask) for masks in per_strategy_masks for mask in masks]
    if len(flattened) != len(set(flattened)):
        raise AssertionError("Current strategies are not composition-deduplicated.")
    return design


def build_source_plan(
    od_readings: dict[str, float | None],
    active_labels: list[str],
    prep_volume_ul: float,
    reading_dilution_factor: float,
    target_working_od: float,
) -> pd.DataFrame:
    rows = []
    for index, label in enumerate(active_labels, 1):
        reading = od_readings.get(label, np.nan)
        reading = np.nan if reading is None else float(reading)
        actual_od = reading * reading_dilution_factor if np.isfinite(reading) else np.nan
        if np.isfinite(actual_od) and actual_od > 0:
            culture_ul = (
                target_working_od * prep_volume_ul / actual_od
                if actual_od > target_working_od
                else prep_volume_ul
            )
            diluent_ul = prep_volume_ul - culture_ul
            working_od = actual_od * culture_ul / prep_volume_ul
            dilution_factor = prep_volume_ul / culture_ul
        else:
            culture_ul = diluent_ul = working_od = dilution_factor = np.nan
        rows.append({
            "lab_marking": label,
            "Source Well": f"A{index}",
            "OD600_reading_at_1_to_10": reading,
            "actual_culture_OD600": actual_od,
            "dilution_factor": dilution_factor,
            "culture_uL_for_working_prep": culture_ul,
            "diluent_uL_for_working_prep": diluent_ul,
            "working_source_OD600": working_od,
            "source_load_uL": prep_volume_ul,
        })
    return pd.DataFrame(rows)


def _transfer_rows(
    design: pd.DataFrame,
    source_plan: pd.DataFrame,
    source_plate_name: str,
    total_biomass_target: float,
    fixed_member_biomass_target: float,
    dosing_mode: str,
    droplet_nl: int,
) -> pd.DataFrame:
    source = source_plan.set_index("lab_marking")
    rows = []
    for destination in design.itertuples(index=False):
        members = destination.members.split("|")
        member_target = (
            total_biomass_target / len(members)
            if dosing_mode == "equal_total_biomass"
            else fixed_member_biomass_target
        )
        for label in members:
            measured_od = float(source.loc[label, "working_source_OD600"])
            calculated_nl = 1000.0 * member_target / measured_od
            executed_nl = int(np.floor(calculated_nl / droplet_nl + 0.5) * droplet_nl)
            rows.append({
                "Source Plate Name": source_plate_name,
                "Source Well": source.loc[label, "Source Well"],
                "Destination Plate Name": getattr(destination, "_2"),
                "Destination Well": getattr(destination, "_3"),
                "Transfer Volume": executed_nl,
                "Source Strategy": destination[0],
                "community_id": destination.community_id,
                "technical_replicate": destination.technical_replicate,
                "community_size": destination.community_size,
                "lab_marking": label,
                "working_source_OD600": measured_od,
                "dosing_mode": dosing_mode,
                "target_member_OD600_uL": member_target,
                "calculated_transfer_nL": calculated_nl,
                "executed_transfer_nL": executed_nl,
                "realized_member_OD600_uL": measured_od * executed_nl / 1000.0,
            })
    return pd.DataFrame(rows)


def make_picklists(
    design: pd.DataFrame,
    source_plan: pd.DataFrame,
    source_plate_name: str,
    max_destination_inoculum_nl: int,
    source_load_ul: float,
    source_dead_volume_ul: float,
    droplet_nl: int,
    full_community_member_biomass: float,
    active_labels: list[str],
) -> tuple[bool, str, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    invalid = source_plan.loc[
        source_plan["working_source_OD600"].isna()
        | source_plan["working_source_OD600"].le(0),
        "lab_marking",
    ].tolist()
    if invalid:
        message = "OD input required for: " + ", ".join(invalid)
        return False, message, pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    total_biomass_target = full_community_member_biomass * len(active_labels)
    transfers = _transfer_rows(
        design, source_plan, source_plate_name, total_biomass_target,
        full_community_member_biomass, "equal_total_biomass", droplet_nl,
    )
    well_audit = transfers.groupby(
        ["Source Strategy", "Destination Plate Name", "Destination Well", "community_id", "technical_replicate", "community_size", "dosing_mode"],
        as_index=False,
    ).agg(
        total_inoculum_nL=("executed_transfer_nL", "sum"),
        realized_total_OD600_uL=("realized_member_OD600_uL", "sum"),
    )
    equal_source_demand = transfers.groupby(
        ["Destination Plate Name", "lab_marking", "Source Well"], as_index=False
    ).agg(transfer_demand_nL=("executed_transfer_nL", "sum"))
    usable_source_ul = source_load_ul - source_dead_volume_ul
    equal_source_demand_exceeds_usable = (
        equal_source_demand["transfer_demand_nL"].max() / 1000.0 > usable_source_ul
    )
    equal_destination_exceeds_cap = (
        well_audit["total_inoculum_nL"].max() > max_destination_inoculum_nl
    )
    if equal_source_demand_exceeds_usable or equal_destination_exceeds_cap:
        transfers = _transfer_rows(
            design, source_plan, source_plate_name, total_biomass_target,
            full_community_member_biomass, "fixed_per_member_biomass", droplet_nl,
        )
        well_audit = transfers.groupby(
            ["Source Strategy", "Destination Plate Name", "Destination Well", "community_id", "technical_replicate", "community_size", "dosing_mode"],
            as_index=False,
        ).agg(
            total_inoculum_nL=("executed_transfer_nL", "sum"),
            realized_total_OD600_uL=("realized_member_OD600_uL", "sum"),
        )

    mode = str(transfers["dosing_mode"].iloc[0])
    well_audit["target_total_OD600_uL"] = (
        total_biomass_target
        if mode == "equal_total_biomass"
        else full_community_member_biomass * well_audit["community_size"]
    )
    well_audit["biomass_error_pct"] = 100.0 * (
        well_audit["realized_total_OD600_uL"] / well_audit["target_total_OD600_uL"] - 1.0
    )

    if (transfers["executed_transfer_nL"] <= 0).any():
        return False, "At least one calculated transfer rounded to zero.", transfers, well_audit, pd.DataFrame()
    if well_audit["total_inoculum_nL"].max() > max_destination_inoculum_nl:
        return False, "Both dosing modes exceed the 15 µL destination inoculum cap.", transfers, well_audit, pd.DataFrame()

    source_demand = transfers.groupby(
        ["Destination Plate Name", "lab_marking", "Source Well"], as_index=False
    ).agg(
        destination_wells=("Destination Well", "nunique"),
        transfer_demand_nL=("executed_transfer_nL", "sum"),
    )
    source_demand["transfer_demand_uL"] = source_demand["transfer_demand_nL"] / 1000.0
    source_demand["source_load_uL"] = source_load_ul
    source_demand["expected_remaining_uL"] = source_load_ul - source_demand["transfer_demand_uL"]
    source_demand["required_dead_volume_uL"] = source_dead_volume_ul
    source_demand["dead_volume_respected"] = (
        source_demand["expected_remaining_uL"] >= source_dead_volume_ul
    )
    if not source_demand["dead_volume_respected"].all():
        return False, "Even fixed-per-member dosing would enter the 15 µL source dead volume.", transfers, well_audit, source_demand

    message = (
        f"Picklists ready in {mode}; maximum destination inoculum "
        f"{well_audit['total_inoculum_nL'].max() / 1000.0:.3f} µL. "
        f"Every 50 µL source fill retains at least {source_dead_volume_ul:.1f} µL. "
        "Refill source row A before each destination plate."
    )
    return True, message, transfers, well_audit, source_demand


def run_adaptive_design(
    *,
    run_label: str,
    od_readings: dict[str, float | None],
    base_dir: Path | str | None = None,
    random_seed: int = 20260902,
    ei_xi: float = 0.002,
    prep_volume_ul: float = 50.0,
    target_working_od: float = 1.20,
    source_dead_volume_ul: float = 15.0,
    max_destination_inoculum_ul: float = 15.0,
    source_plate_name: str = "396",
    excluded_members: set[str] | None = None,
    export: bool = True,
) -> dict[str, object]:
    """Load round 1, fit the GP, choose four plates, and optionally export files."""
    base_dir = Path(base_dir) if base_dir is not None else Path(__file__).resolve().parent
    excluded_members = set() if excluded_members is None else set(excluded_members)
    unknown_exclusions = excluded_members - set(MEMBER_LABELS)
    if unknown_exclusions:
        raise ValueError(f"Unknown excluded members: {sorted(unknown_exclusions)}")
    active_labels = [label for label in MEMBER_LABELS if label not in excluded_members]
    if not active_labels:
        raise ValueError("At least one active member is required.")
    previous_design_path = base_dir / "experiment_setup_output" / "community_design_20260828_13strains.csv"
    readout_path = base_dir / "readout" / "Karl_20260830_201548_MAGNUS_ABS-Lactate-assay-rerun.csv"
    output_dir = base_dir / "experiment_setup_output" / f"{run_label}_adaptive_{len(active_labels)}strains"

    previous_design = pd.read_csv(previous_design_path)
    training, blank_median = load_training_data(previous_design_path, readout_path)
    x_train = training[[f"member_{label}" for label in MEMBER_LABELS]].to_numpy(float)
    y_train = training["response_blank_corrected"].to_numpy(float)
    gp_fit, diagnostics = fit_gp(x_train, y_train)
    diagnostics.update({
        "active_species_count": len(active_labels),
        "excluded_species": "|".join(sorted(excluded_members)) or "none",
        "edge_blank_median": blank_median,
        "training_response": "triplicate median absorbance minus edge-blank median",
        "readout_file": readout_path.name,
        "previous_design_file": previous_design_path.name,
    })
    model_summary = pd.DataFrame([diagnostics])

    predictions = build_candidate_predictions(
        training, gp_fit, ei_xi=ei_xi, active_labels=active_labels
    )
    batches, predictions = select_strategies(
        predictions, previous_design, n_per_strategy=102, random_seed=random_seed,
        active_labels=active_labels,
    )
    design = build_plate_designs(
        batches, previous_design, run_label, random_seed, technical_replicates=3
    )

    strategy_summary = design.groupby(
        ["Source Strategy", "Destination Plate Name"], as_index=False
    ).agg(
        unique_communities=("community_mask", "nunique"),
        destination_wells=("Destination Well", "nunique"),
        min_community_size=("community_size", "min"),
        median_community_size=("community_size", "median"),
        max_community_size=("community_size", "max"),
        exact_repeat_wells=("repeat_status", lambda values: int(values.eq("exact_repeat").sum())),
        maxdiv_backfill_wells=("repeat_status", lambda values: int(values.eq("active_panel_maxdiv_backfill").sum())),
    )

    source_plan = build_source_plan(
        od_readings, active_labels=active_labels, prep_volume_ul=prep_volume_ul,
        reading_dilution_factor=10.0, target_working_od=target_working_od,
    )
    picklist_ready, picklist_message, transfers, well_audit, source_demand = make_picklists(
        design,
        source_plan,
        source_plate_name=source_plate_name,
        max_destination_inoculum_nl=int(max_destination_inoculum_ul * 1000),
        source_load_ul=prep_volume_ul,
        source_dead_volume_ul=source_dead_volume_ul,
        droplet_nl=25,
        full_community_member_biomass=0.100,
        active_labels=active_labels,
    )

    export_paths: list[Path] = []
    if export:
        output_dir.mkdir(parents=True, exist_ok=True)
        frames = {
            f"training_data_{run_label}.csv": training,
            f"gp_model_summary_{run_label}.csv": model_summary,
            f"candidate_predictions_{run_label}.csv": predictions,
            f"sampling_strategy_summary_{run_label}.csv": strategy_summary,
            f"community_design_all_plates_{run_label}.csv": design,
            f"source_plate_map_and_dilution_{run_label}.csv": source_plan,
        }
        for filename, frame in frames.items():
            path = output_dir / filename
            frame.to_csv(path, index=False)
            export_paths.append(path)

        for plate_name, plate_design in design.groupby("Destination Plate Name"):
            safe_name = plate_name.replace("community_", "")
            path = output_dir / f"community_design_{safe_name}.csv"
            plate_design.to_csv(path, index=False)
            export_paths.append(path)

        if picklist_ready:
            transfer_audit_path = output_dir / f"transfer_audit_all_plates_{run_label}.csv"
            transfers.to_csv(transfer_audit_path, index=False)
            export_paths.append(transfer_audit_path)
            well_audit_path = output_dir / f"well_audit_all_plates_{run_label}.csv"
            well_audit.to_csv(well_audit_path, index=False)
            export_paths.append(well_audit_path)
            source_demand_path = output_dir / f"source_demand_by_plate_{run_label}.csv"
            source_demand.to_csv(source_demand_path, index=False)
            export_paths.append(source_demand_path)
            echo_columns = [
                "Source Plate Name", "Source Well", "Destination Plate Name",
                "Destination Well", "Transfer Volume",
            ]
            for plate_name, plate_transfers in transfers.groupby("Destination Plate Name"):
                safe_name = plate_name.replace("community_", "")
                path = output_dir / f"echo_picklist_{safe_name}.csv"
                plate_transfers[echo_columns].to_csv(path, index=False)
                export_paths.append(path)

        manifest = pd.DataFrame({
            "run_label": run_label,
            "artifact": [path.name for path in export_paths],
            "picklist_ready": picklist_ready,
            "picklist_status": picklist_message,
        })
        manifest_path = output_dir / f"run_manifest_{run_label}.csv"
        manifest.to_csv(manifest_path, index=False)
        export_paths.append(manifest_path)

    return {
        "training_data": training,
        "model_summary": model_summary,
        "candidate_predictions": predictions,
        "strategy_summary": strategy_summary,
        "design": design,
        "source_plan": source_plan,
        "picklist_ready": picklist_ready,
        "picklist_message": picklist_message,
        "transfer_audit": transfers,
        "well_audit": well_audit,
        "source_demand": source_demand,
        "output_dir": output_dir,
        "export_paths": export_paths,
    }


if __name__ == "__main__":
    current_od_readings = {
        "256": 0.53, "75": 0.61, "269": 0.06, "294": 0.28,
        "275": 1.22, "41": 0.58, "158": 0.81, "145": 0.53,
        "200": 0.22, "196": 0.08, "291": 0.13, "296": 0.00,
        "264": 0.76,
    }
    result = run_adaptive_design(
        run_label="20260902",
        od_readings=current_od_readings,
        excluded_members={"296"},
        export=True,
    )
    print(result["model_summary"].to_string(index=False))
    print(result["strategy_summary"].to_string(index=False))
    print(result["picklist_message"])
    print(f"Exported design/model artifacts to {result['output_dir']}")


# ---------------- Appended 2026-09-09: independent acquisition paths ----------------
# The 2026-09-02 definitions and direct-script entry point above are historical.
# Execute the new dated notebook section or acquisition_round_20260909.py for run 3.
def run_acquisition_round_20260909(**kwargs):
    """Continue EI/TM models and random/max-diversity paths without rewriting run 2."""
    from acquisition_round_20260909 import run_round
    return run_round(**kwargs)


# ---------------- Appended 2026-09-16: corrected history, four source rows ----------------
# Edit ONLY today's values here. These are OD600 readings at 1:10 dilution;
# the new round multiplies by 10. None means pending, not zero growth.
# Source columns 1–12 repeat in rows A/B/C/D for destination plates 1/2/3/4.
# 50 µL per source well; 200 µL of working suspension per species across four rows.
# 296 stays excluded from today's fixed 12-species panel.
# Added validation reservation: plate 2 has 94 new max-diversity communities plus
# 3 community repeats and their 5 constituent monocultures, all in triplicate.
# The original four-member candidate lost 196 in one well; repeat both the full
# community and its recorded three-member version, plus the intact 294+275 pair.
OD600_READING_AT_1_TO_10_20260916 = {
    "256": 0.62, "75": 0.64, "269": 0.07, "294": 0.41,
    "275": 1.36, "41": 0.67, "158": 0.60, "145": 0.74,
    "200": 0.27, "196": 0.09, "291": 0.15, "264": 0.59,
}


def run_acquisition_round_20260916(od_readings=None, **kwargs):
    """Fit transfer-corrected paths; missing ODs export a design but no Echo files.

    Run acquisition_round_20260916.py directly for this round. The earlier
    adaptive_run.py direct-script entry point remains the historical Sep 2 run.
    """
    from acquisition_round_20260916 import run_round
    readings = OD600_READING_AT_1_TO_10_20260916 if od_readings is None else od_readings
    return run_round(od_readings=readings, **kwargs)


# ---------------- Appended 2026-09-23: one-plate analytical diagnostic round ----------------
# No acquisition, model fitting, or changes to previous ODs/results.
# Enter TODAY'S OD600 readings at 1:10 dilution. None = not measured yet.
# Source plate 396, A1-A12, 50 uL EACH; no need for four source rows this round.
OD600_READING_AT_1_TO_10_20260923 = {
    "256": None, "75": None, "269": None, "294": None,
    "275": None, "41": None, "158": None, "145": None,
    "200": None, "196": None, "291": None, "264": None,
}


def run_analytical_round_20260923(od_readings=None, **kwargs):
    """One culture plate and one later assay plate; no Echo export until ODs exist.

    Run analytical_round_20260923.py or the final notebook cell for this round.
    Running adaptive_run.py directly still executes its historical Sep 2 entry point.
    """
    from analytical_round_20260923 import run_round
    readings = OD600_READING_AT_1_TO_10_20260923 if od_readings is None else od_readings
    return run_round(od_readings=readings, **kwargs)
