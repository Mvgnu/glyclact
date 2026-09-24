#!/usr/bin/env python3
"""Assign structural fold to every sequence by HMM, for the panel and the candidates.

Why: the first pass of the sequence screen produced 186 genomes that gained an L-LDH
relative to the annotation, and the cause was a mis-folded *reference*. Swiss-Prot
Q07251 (Cupriavidus necator, gene `ldh`) carries EC 1.1.1.27 but its only Pfam domain is
PF02615 `Ldh_2` -- the malate/L-lactate dehydrogenase family, a different fold from the
canonical L-LDH. It recruited 141 Enterobacteriaceae `ybiC` proteins as false L-LDH. So
the reference set built to catch annotation artifacts had one inside it.

A label cannot be trusted on either side of the comparison. Fold can, because it is
measured rather than transferred. Three architectures matter here and they are mutually
exclusive:

  Ldh_1     PF00056 + PF02866    NAD-dependent L-LDH *and* malate dehydrogenase
  2-Hacid   PF00389 + PF02826    D-LDH and its whole decoy family (serA, ghrB, pdxB, ...)
  Ldh_2     PF02615              LDH2/MDH2 family; not the fermentative L-LDH

So a D-LDH call requires the 2-Hacid architecture, an L-LDH call requires Ldh_1, and
anything scoring Ldh_2 is excluded from both regardless of what its best hit was named.
Within Ldh_1, L-LDH and MDH share the fold and are separated by best-hit scoring against
a panel that holds both -- fold narrows the question rather than answering it.

Writes references/ldh_panel/panel_folds.tsv and, given a query FASTA, a folds table for
the candidates.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PANEL_DIR = ROOT / "references" / "ldh_panel"

HMMS = {
    "Ldh_1_N": PANEL_DIR / "PF00056.hmm",
    "Ldh_1_C": PANEL_DIR / "PF02866.hmm",
    "Ldh_2": PANEL_DIR / "PF02615.hmm",
    "2-Hacid_dh": PANEL_DIR / "PF00389.hmm",
    "2-Hacid_dh_C": PANEL_DIR / "PF02826.hmm",
}

# Fold is decided by which architecture scores best, not by presence alone: the N- and
# C-terminal halves of Ldh_1 cross-hit 2-Hacid weakly in some proteins.
ARCHITECTURE = {
    "Ldh_1": ("Ldh_1_N", "Ldh_1_C"),
    "2-Hacid": ("2-Hacid_dh", "2-Hacid_dh_C"),
    "Ldh_2": ("Ldh_2",),
}

DOMAIN_EVALUE = 1e-5


def hmmsearch(hmm: Path, fasta: Path, threads: int) -> dict[str, float]:
    """Best per-sequence bitscore for one HMM."""
    proc = subprocess.run(
        ["hmmsearch", "--cpu", str(threads), "-E", str(DOMAIN_EVALUE),
         "--tblout", "/dev/stdout", "-o", "/dev/null", str(hmm), str(fasta)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"hmmsearch failed on {hmm.name}:\n{proc.stderr[-1500:]}")
    best: dict[str, float] = {}
    for line in proc.stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        name, score = parts[0], float(parts[5])
        if score > best.get(name, -1e9):
            best[name] = score
    return best


def assign(fasta: Path, threads: int) -> dict[str, tuple[str, dict[str, float]]]:
    scores: dict[str, dict[str, float]] = defaultdict(dict)
    for label, hmm in HMMS.items():
        for name, score in hmmsearch(hmm, fasta, threads).items():
            scores[name][label] = score

    result = {}
    for name, s in scores.items():
        arch_score = {}
        for arch, doms in ARCHITECTURE.items():
            # Sum the constituent domains; a two-domain architecture needs both to be
            # present or it loses to the single-domain one on a partial match.
            if all(d in s for d in doms):
                arch_score[arch] = sum(s[d] for d in doms)
            elif any(d in s for d in doms):
                arch_score[arch] = max(s[d] for d in doms if d in s) * 0.5
        if not arch_score:
            continue
        fold = max(arch_score, key=arch_score.get)
        result[name] = (fold, s)
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", type=Path,
                    help="candidate FASTA; if omitted only the panel is profiled")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    targets = [(PANEL_DIR / "ldh_panel.faa", PANEL_DIR / "panel_folds.tsv")]
    if args.query:
        targets.append((args.query, args.out or args.query.with_suffix(".folds.tsv")))

    for fasta, out in targets:
        folds = assign(fasta, args.threads)
        with out.open("w") as fh:
            fh.write("sequence\tfold\t" + "\t".join(HMMS) + "\n")
            for name, (fold, s) in sorted(folds.items()):
                fh.write(name + "\t" + fold + "\t"
                         + "\t".join(f"{s.get(k, 0):.1f}" for k in HMMS) + "\n")
        print(f"wrote {out} ({len(folds)} sequences)", file=sys.stderr)

        if fasta.name == "ldh_panel.faa":
            # A reference whose fold contradicts its class is the failure mode that
            # caused this script to exist, so it is reported rather than silently fixed.
            expected = {
                "L_LDH": "Ldh_1", "MDH": "Ldh_1",
                "D_LDH": "2-Hacid", "SERA": "2-Hacid", "GHR": "2-Hacid",
                "PDXB": "2-Hacid", "D2HA": "2-Hacid", "FDH_NAD": "2-Hacid",
            }
            bad = defaultdict(list)
            for name, (fold, _) in folds.items():
                cls = name.split("|")[0]
                if cls in expected and fold != expected[cls]:
                    bad[cls].append((name, fold))
            print("\npanel references whose measured fold contradicts their EC class:",
                  file=sys.stderr)
            if not bad:
                print("  none", file=sys.stderr)
            for cls, items in sorted(bad.items()):
                print(f"  {cls}: {len(items)}", file=sys.stderr)
                for name, fold in items[:8]:
                    print(f"      {name}  -> {fold}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
