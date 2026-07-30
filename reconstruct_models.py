#!/usr/bin/env python3
"""Reconstruct the local candidate panel with CarveMe from Prokka proteins."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from pipeline.gem_reconstruction.reconstruct import patch_carveme_pandas_compatibility


def reconstruct(
    manifest_path: Path,
    output_dir: Path,
    gapfill_medium: str,
    selected_genomes: set[str] | None = None,
) -> None:
    from carveme.cli.carve import maincall
    from reframed import set_default_solver

    patch_carveme_pandas_compatibility()
    set_default_solver("scip")
    output_dir.mkdir(parents=True, exist_ok=True)
    with manifest_path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if selected_genomes:
        rows = [row for row in rows if row["genome"] in selected_genomes]
        missing = selected_genomes - {row["genome"] for row in rows}
        if missing:
            raise ValueError(f"Genome(s) absent from manifest: {sorted(missing)}")
    for index, row in enumerate(rows, 1):
        genome = row["genome"]
        destination = output_dir / f"{genome}.xml"
        if destination.exists() and destination.stat().st_size > 0:
            print(f"[{index}/{len(rows)}] {genome}: existing model retained", flush=True)
            continue
        # A genome that has been deliberately excluded still sits in the manifest, and its
        # model lives in models/excluded/. Without this check a later run rebuilds it
        # straight back into models/ -- which silently undid seven exclusions once.
        if (output_dir / "excluded" / f"{genome}.xml").exists():
            print(f"[{index}/{len(rows)}] {genome}: excluded model exists, not rebuilding",
                  flush=True)
            continue
        print(f"[{index}/{len(rows)}] {genome}: reconstructing", flush=True)
        maincall(
            inputfile=row["protein_fasta"],
            input_type="protein",
            outputfile=str(destination),
            universe="bacteria",
            flavor="fbc2",
            gapfill=gapfill_medium or None,
            default_score=-1.0,
            uptake_score=0.0,
            soft_score=1.0,
            ref_score=0.0,
            verbose=True,
        )
        if not destination.exists():
            raise RuntimeError(f"CarveMe failed to write {destination}")


def parse_args() -> argparse.Namespace:
    root = Path(__file__).parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=root / "candidate_genomes" / "manifest.tsv")
    parser.add_argument("--output-dir", type=Path, default=root / "models")
    # M9 is glucose-based. It repairs generic biomass feasibility without forcing the target
    # glycerol or DHA phenotype into a model during gap filling.
    parser.add_argument("--gapfill-medium", default="M9")
    parser.add_argument("--genome", action="append", help="Reconstruct only this genome (repeatable).")
    return parser.parse_args()


if __name__ == "__main__":
    options = parse_args()
    reconstruct(
        options.manifest,
        options.output_dir,
        options.gapfill_medium,
        set(options.genome) if options.genome else None,
    )
