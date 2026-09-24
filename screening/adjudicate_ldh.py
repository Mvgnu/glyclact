#!/usr/bin/env python3
"""Decide L-LDH / D-LDH membership by sequence, not by annotation transfer.

Why this exists
---------------
Pool v4 rests on per-genome L-LDH and D-LDH counts, and two tables in this repo
disagree about those counts in 19 of 260 genomes. Worse, both tables get their answer
the same way -- from eggNOG's EC/KO/name fields -- so agreement between them would not
have been evidence either. A D-LDH call from an automatic annotator is a *family*
assignment: `ldhD` sits in the D-isomer-specific 2-hydroxyacid dehydrogenase family
alongside `serA`, `ghrA`/`ghrB`, `pdxB`, glycerate dehydrogenase, NAD-formate
dehydrogenase and D-2-hydroxyisocaproate dehydrogenase, all of which share its fold,
its PFAM architecture (PF00389 + PF02826) and its His/Glu/Arg catalytic triad. Several
never see pyruvate in vivo.

That matters here because the purity gate is absolute: one real D-LDH takes a batch
off-spec, but one *spurious* D-LDH call silently deletes a strain the pool needs.

What this does
--------------
Four lines of evidence per candidate protein, each recorded separately so a
call can be inspected rather than trusted:

0. **Fold, measured by HMM, on both sides of the comparison.** This gate exists because
   the first version of this screen produced 186 genomes that gained an L-LDH, and the
   cause was a bad *reference*: Swiss-Prot Q07251 (Cupriavidus necator `ldh`) carries
   EC 1.1.1.27 but its only Pfam domain is PF02615 `Ldh_2`, a different fold from the
   canonical L-LDH, and it recruited 141 Enterobacteriaceae `ybiC` proteins. A reference
   panel assembled from EC labels can therefore carry the same disease it was built to
   diagnose. Every panel member and every candidate is now fold-typed by hmmsearch
   (see assign_folds.py); a reference whose measured fold contradicts its EC class is
   demoted to a decoy, and a candidate must carry the right architecture -- Ldh_1 for
   L-LDH, 2-Hacid for D-LDH -- before either call is available to it.

1. **Best-hit class.** Every protein is searched against a Swiss-Prot panel holding the
   targets *and* their decoys together (see build_ldh_references.py). The class of the
   best hit names the protein. A protein that prefers SerA to any characterized LdhD is
   a family member, and this is the single strongest discriminator.

2. **Class margin.** Bitscore of the best target hit minus the best decoy hit. A real
   ortholog sits well clear of the decoys; a family member sits between them or on the
   wrong side. Reported as a number, not a verdict.

3. **Reciprocal best hit.** The reference that a candidate matched is searched back
   against that genome's own proteome. RBH holds only if the candidate is the top
   returning protein. This is the standard orthology test and the one that separates
   "closest member of the family I have" from "ortholog".

Writes results/screening/ldh_adjudication.tsv (per protein) and
results/screening/ldh_calls.tsv (per genome).
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PANEL_DIR = ROOT / "references" / "ldh_panel"
PANEL_FAA = PANEL_DIR / "ldh_panel.faa"
PANEL_DMND = PANEL_DIR / "ldh_panel.dmnd"
SCORES = ROOT / "results" / "screening" / "genome_scores.tsv"
OUT_DIR = ROOT / "results" / "screening"

TARGETS = {"L_LDH", "D_LDH"}
DECOYS = {"MDH", "SERA", "GHR", "PDXB", "D2HA", "FDH_NAD", "FOLD_MISMATCH"}
OTHER = {"LLDD", "DLD", "LOX", "LARA", "LAR_MAT"}
REPORTABLE = TARGETS | OTHER

# The architecture each class must actually have. Classes absent here are not
# fold-checked: the lactate consumers and the lar system are unrelated folds whose
# identification does not depend on this axis.
EXPECTED_FOLD = {
    "L_LDH": "Ldh_1", "MDH": "Ldh_1",
    "D_LDH": "2-Hacid", "SERA": "2-Hacid", "GHR": "2-Hacid",
    "PDXB": "2-Hacid", "D2HA": "2-Hacid", "FDH_NAD": "2-Hacid",
}

# A hit weaker than this is not evidence of anything.
MIN_IDENTITY = 30.0
MIN_COVERAGE = 0.60
MAX_EVALUE = 1e-20

FIELDS = ("qseqid sseqid pident length qlen slen evalue bitscore "
          "qstart qend sstart send")


def sh(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if proc.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd[:3])} failed:\n{proc.stderr[-2000:]}")
    return proc


def load_genomes() -> list[dict]:
    with SCORES.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build_combined_proteome(genomes: list[dict], dest: Path) -> dict[str, str]:
    """One FASTA for all genomes, IDs prefixed so hits stay attributable."""
    protein_to_genome: dict[str, str] = {}
    with dest.open("w") as out:
        for row in genomes:
            gid, faa = row["genome"], row["protein_fasta"]
            if not faa or not Path(faa).exists():
                print(f"  missing proteome: {gid}", file=sys.stderr)
                continue
            with open(faa) as fh:
                for line in fh:
                    if line.startswith(">"):
                        pid = line[1:].split()[0]
                        tagged = f"{gid}~{pid}"
                        protein_to_genome[tagged] = gid
                        out.write(f">{tagged}\n")
                    else:
                        out.write(line)
    return protein_to_genome


def diamond_search(query: Path, db: Path, out: Path, k: int, threads: int) -> None:
    sh([
        "diamond", "blastp",
        "-q", str(query), "-d", str(db), "-o", str(out),
        "--outfmt", "6", *FIELDS.split(),
        "--max-target-seqs", str(k),
        "--evalue", "1e-5",
        "--very-sensitive",
        "--threads", str(threads),
        "--quiet",
    ])


def read_hits(path: Path):
    with path.open() as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            yield {
                "q": p[0], "s": p[1],
                "pident": float(p[2]), "length": int(p[3]),
                "qlen": int(p[4]), "slen": int(p[5]),
                "evalue": float(p[6]), "bits": float(p[7]),
                "qstart": int(p[8]), "qend": int(p[9]),
                "sstart": int(p[10]), "send": int(p[11]),
            }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=os.environ.get("LDH_WORK",
                                                    str(ROOT / "work" / "ldh")),
                    help="scratch for the concatenated proteome and search tables; "
                         "large, and gitignored")
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--reuse", action="store_true",
                    help="reuse existing diamond output if present")
    args = ap.parse_args()

    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    genomes = load_genomes()
    print(f"genomes: {len(genomes)}", file=sys.stderr)

    combined = work / "all_proteins.faa"
    if not (args.reuse and combined.exists()):
        print("concatenating proteomes ...", file=sys.stderr)
        build_combined_proteome(genomes, combined)
    n_prot = sum(1 for line in combined.open() if line.startswith(">"))
    print(f"proteins: {n_prot}", file=sys.stderr)

    # Forward: every genome protein against the panel.
    fwd = work / "forward.tsv"
    if not (args.reuse and fwd.exists()):
        print("forward search (proteome -> panel) ...", file=sys.stderr)
        diamond_search(combined, PANEL_DMND, fwd, k=60, threads=args.threads)

    # Reverse: the panel back against every genome protein, for the RBH test.
    all_db = work / "all_proteins"
    rev = work / "reverse.tsv"
    if not (args.reuse and rev.exists()):
        print("building proteome database ...", file=sys.stderr)
        sh(["diamond", "makedb", "--in", str(combined), "-d", str(all_db),
            "--threads", str(args.threads), "--quiet"])
        print("reverse search (panel -> proteome) ...", file=sys.stderr)
        diamond_search(PANEL_FAA, all_db.with_suffix(".dmnd"), rev,
                       k=25000, threads=args.threads)

    # --- fold-type the panel and the candidates ------------------------------
    import assign_folds

    panel_folds_path = PANEL_DIR / "panel_folds.tsv"
    if not panel_folds_path.exists():
        raise SystemExit("run assign_folds.py first to profile the panel")
    panel_fold = {}
    with panel_folds_path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            panel_fold[row["sequence"]] = row["fold"]

    # Demote references whose measured fold contradicts their EC class.
    ref_effective: dict[str, str] = {}
    demoted = 0
    for name, fold in panel_fold.items():
        cls = name.split("|")[0]
        if cls in EXPECTED_FOLD and fold != EXPECTED_FOLD[cls]:
            ref_effective[name] = "FOLD_MISMATCH"
            demoted += 1
        else:
            ref_effective[name] = cls
    print(f"panel references demoted for fold mismatch: {demoted}", file=sys.stderr)

    # Only proteins that hit the panel need fold-typing.
    hit_ids = {line.split("\t", 1)[0] for line in fwd.open()}
    cand_faa = work / "candidates.faa"
    if not (args.reuse and cand_faa.exists()):
        keep = False
        with cand_faa.open("w") as out, combined.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    keep = line[1:].strip().split()[0] in hit_ids
                if keep:
                    out.write(line)
    print(f"fold-typing {len(hit_ids)} candidate proteins ...", file=sys.stderr)
    cand_fold = {name: fold
                 for name, (fold, _) in assign_folds.assign(cand_faa, args.threads).items()}

    # --- collate forward hits -------------------------------------------------
    # per protein: best bitscore per class, and the best hit overall
    per_protein: dict[str, dict] = {}
    for h in read_hits(fwd):
        cls, acc, gene, org = h["s"].split("|")
        cls = ref_effective.get(h["s"], cls)
        cov = h["length"] / max(h["qlen"], h["slen"])
        if h["pident"] < MIN_IDENTITY or cov < MIN_COVERAGE or h["evalue"] > MAX_EVALUE:
            continue
        rec = per_protein.setdefault(h["q"], {"best": None, "by_class": {}})
        cur = rec["by_class"].get(cls)
        if cur is None or h["bits"] > cur["bits"]:
            rec["by_class"][cls] = {
                "bits": h["bits"], "acc": acc, "gene": gene, "org": org,
                "pident": h["pident"], "cov": cov, "evalue": h["evalue"],
            }
        if rec["best"] is None or h["bits"] > rec["best"]["bits"]:
            rec["best"] = {
                "cls": cls, "bits": h["bits"], "acc": acc, "gene": gene,
                "org": org, "pident": h["pident"], "cov": cov,
                "evalue": h["evalue"],
            }

    print(f"proteins with a panel hit: {len(per_protein)}", file=sys.stderr)

    # --- collate reverse hits for RBH ----------------------------------------
    # For each (reference accession, genome): which protein scores highest?
    rbh_top: dict[tuple[str, str], tuple[str, float]] = {}
    for h in read_hits(rev):
        cls, acc, gene, org = h["q"].split("|")
        gid = h["s"].split("~", 1)[0]
        key = (acc, gid)
        cur = rbh_top.get(key)
        if cur is None or h["bits"] > cur[1]:
            rbh_top[key] = (h["s"], h["bits"])

    # --- per-protein adjudication --------------------------------------------
    rows = []
    for pid, rec in per_protein.items():
        best = rec["best"]
        if best["cls"] not in REPORTABLE:
            continue
        gid = pid.split("~", 1)[0]
        by_class = rec["by_class"]

        best_decoy = max(
            (v["bits"] for k, v in by_class.items() if k in DECOYS), default=0.0
        )
        best_decoy_cls = max(
            ((v["bits"], k) for k, v in by_class.items() if k in DECOYS),
            default=(0.0, "none"),
        )[1]
        margin = best["bits"] - best_decoy

        # The other isomer's best score, so an L call can be checked against D directly.
        other = "D_LDH" if best["cls"] == "L_LDH" else "L_LDH"
        other_bits = by_class.get(other, {}).get("bits", 0.0)

        rbh_protein, _ = rbh_top.get((best["acc"], gid), ("", 0.0))
        rbh = rbh_protein == pid

        # Fold gate, query side. An isomer call is only available to a protein that
        # actually has the architecture; without it the call is recorded but marked,
        # never counted.
        fold = cand_fold.get(pid, "none")
        want = EXPECTED_FOLD.get(best["cls"])
        fold_ok = (want is None) or (fold == want)

        rows.append({
            "genome": gid,
            "protein": pid.split("~", 1)[1],
            "call": best["cls"],
            "fold": fold,
            "fold_expected": want or "n/a",
            "fold_ok": "1" if fold_ok else "0",
            "best_ref": best["acc"],
            "best_ref_gene": best["gene"],
            "best_ref_organism": best["org"],
            "pident": f"{best['pident']:.1f}",
            "coverage": f"{best['cov']:.2f}",
            "evalue": f"{best['evalue']:.1e}",
            "bitscore": f"{best['bits']:.1f}",
            "best_decoy_class": best_decoy_cls,
            "best_decoy_bits": f"{best_decoy:.1f}",
            "decoy_margin": f"{margin:.1f}",
            "other_isomer_bits": f"{other_bits:.1f}",
            "reciprocal_best_hit": "1" if rbh else "0",
            "rbh_winner": "" if rbh else rbh_protein.split("~", 1)[-1],
        })

    rows.sort(key=lambda r: (r["genome"], r["call"], -float(r["bitscore"])))
    per_protein_path = OUT_DIR / "ldh_adjudication.tsv"
    with per_protein_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {per_protein_path} ({len(rows)} proteins)", file=sys.stderr)

    # --- per-genome roll-up ---------------------------------------------------
    # A call counts as confirmed only if it is both best-hit and reciprocal-best-hit.
    by_genome: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_genome[r["genome"]][r["call"]].append(r)

    meta = {row["genome"]: row for row in genomes}
    out_rows = []
    for gid in sorted(meta):
        g = by_genome.get(gid, {})

        def counts(cls):
            hits = [h for h in g.get(cls, []) if h["fold_ok"] == "1"]
            confirmed = [h for h in hits if h["reciprocal_best_hit"] == "1"]
            return len(hits), len(confirmed)

        l_all, l_conf = counts("L_LDH")
        d_all, d_conf = counts("D_LDH")
        lldd_all, lldd_conf = counts("LLDD")
        dld_all, dld_conf = counts("DLD")
        lox_all, _ = counts("LOX")
        lara_all, lara_conf = counts("LARA")
        larmat_all, _ = counts("LAR_MAT")

        d_details = "; ".join(
            f"{h['protein']}:{h['best_ref_gene']}({h['pident']}%,"
            f"m{h['decoy_margin']},rbh{h['reciprocal_best_hit']})"
            for h in g.get("D_LDH", [])
        )
        l_details = "; ".join(
            f"{h['protein']}:{h['best_ref_gene']}({h['pident']}%,"
            f"m{h['decoy_margin']},rbh{h['reciprocal_best_hit']})"
            for h in g.get("L_LDH", [])
        )

        out_rows.append({
            "genome": gid,
            "genus": meta[gid]["genus"],
            "l_ldh_seq": l_conf,
            "d_ldh_seq": d_conf,
            "l_ldh_besthit_only": l_all,
            "d_ldh_besthit_only": d_all,
            "l_ldh_annotation": meta[gid]["l_ldh_count"],
            "d_ldh_annotation": meta[gid]["d_ldh_count"],
            "agrees_with_annotation":
                "1" if (str(l_conf) == meta[gid]["l_ldh_count"]
                        and str(d_conf) == meta[gid]["d_ldh_count"]) else "0",
            "lldD_eats_L": lldd_conf,
            "dld_eats_D": dld_conf,
            "lactate_oxidase": lox_all,
            "larA": lara_conf,
            "lar_maturation": larmat_all,
            "d_ldh_evidence": d_details,
            "l_ldh_evidence": l_details,
        })

    calls_path = OUT_DIR / "ldh_calls.tsv"
    with calls_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {calls_path} ({len(out_rows)} genomes)", file=sys.stderr)

    disagree = sum(1 for r in out_rows if r["agrees_with_annotation"] == "0")
    print(f"\ngenomes where sequence call differs from annotation: "
          f"{disagree}/{len(out_rows)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
