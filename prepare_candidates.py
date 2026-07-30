#!/usr/bin/env python3
"""Copy a declared candidate panel into stable, standardized local inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_one(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one {pattern!r} below {directory}, found {len(matches)}: {matches}"
        )
    return matches[0]


def prepare(config_path: Path, output_dir: Path, manifest_path: Path) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    source_root = Path(config["source_root"])
    if not source_root.is_dir():
        raise FileNotFoundError(
            f"Candidate source root is unavailable: {source_root}. Mount the strain-library drive."
        )

    rows: list[dict[str, str]] = []
    output_dir.mkdir(parents=True, exist_ok=True)
    for entry in config["candidates"]:
        genome = str(entry["genome"])
        source = source_root / genome
        if not source.is_dir():
            raise FileNotFoundError(f"Missing candidate folder: {source}")
        nucleotide = resolve_one(source, "prokka/*.fna")
        protein = resolve_one(source, "prokka/*.faa")
        annotation = resolve_one(source, "eggnog.emapper.annotations")
        destination = output_dir / genome
        destination.mkdir(parents=True, exist_ok=True)

        copied = {
            "nucleotide_fasta": (nucleotide, destination / f"{genome}.fna"),
            "protein_fasta": (protein, destination / f"{genome}.faa"),
            "eggnog_annotations": (annotation, destination / "eggnog.emapper.annotations"),
        }
        for _, (source_path, destination_path) in copied.items():
            shutil.copy2(source_path, destination_path)

        for taxonomy_name in ("unique_taxa.txt", "unique_taxa_percentage.txt"):
            taxonomy_path = source / "diamond_megan" / taxonomy_name
            if taxonomy_path.exists():
                shutil.copy2(taxonomy_path, destination / taxonomy_name)

        # The screen reads assembly QC from quast_output/report.tsv. Without it the
        # local copies report no_quast and lose the independent size/GC check on the
        # taxonomy call.
        quast_report = source / "quast_output" / "report.tsv"
        if quast_report.exists():
            (destination / "quast_output").mkdir(exist_ok=True)
            shutil.copy2(quast_report, destination / "quast_output" / "report.tsv")

        row = {
            "genome": genome,
            "role": str(entry["role"]),
            "rationale": str(entry["rationale"]),
            "source_folder": str(source),
        }
        for label, (_, destination_path) in copied.items():
            row[label] = str(destination_path)
            row[f"{label}_sha256"] = sha256(destination_path)
        rows.append(row)

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("candidate_panel.yaml"))
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).with_name("candidate_genomes")
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(__file__).with_name("candidate_genomes") / "manifest.tsv",
    )
    return parser.parse_args()


if __name__ == "__main__":
    options = parse_args()
    prepare(options.config, options.output_dir, options.manifest)
