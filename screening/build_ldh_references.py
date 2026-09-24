#!/usr/bin/env python3
"""Fetch a curated reference set for adjudicating lactate dehydrogenase calls.

Annotation-level LDH calls are not trustworthy for this project, and the reason is
structural rather than accidental.

**D-LDH sits inside a large family of near neighbours.** `ldhD` (EC 1.1.1.28) belongs to
the D-isomer-specific 2-hydroxyacid dehydrogenase family (PF00389 + PF02826). So do
D-3-phosphoglycerate dehydrogenase (`serA`), glyoxylate and hydroxypyruvate reductases
(`ghrA`/`ghrB`), glycerate dehydrogenase, erythronate-4-phosphate dehydrogenase (`pdxB`),
NAD-dependent formate dehydrogenase, D-mandelate dehydrogenase and
D-2-hydroxyisocaproate dehydrogenase. These share fold, PFAM architecture and the
catalytic His/Glu/Arg triad; several never touch pyruvate in vivo. An automatic
annotator handed a family member will frequently emit "D-lactate dehydrogenase", which
means a D-LDH call from eggNOG is a family assignment, not an orthology assignment.

**L-LDH has one dominant near neighbour**, malate dehydrogenase (EC 1.1.1.37), same fold
and same PFAM architecture (PF00056 + PF02866). The existing screen already handles this
by EC/KO veto; sequence-level best-hit assignment handles it without needing the EC field
to be present.

The fix in both directions is the same: compare against *characterized* sequences and let
the best hit name the class. This script assembles that reference set from Swiss-Prot,
which is manually curated, and restricts to bacteria so the identity thresholds mean
something.

Writes references/ldh_panel/ .
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "references" / "ldh_panel"
STREAM = "https://rest.uniprot.org/uniprotkb/stream"

BACTERIA = "taxonomy_id:2"

# Each class is (label, query, why it is here).
#
# TARGET  -- what a producer must have
# DECOY   -- what gets miscalled as the target
# CONSUME -- eats lactate; matters for titre and for optical purity via lldD
# PURITY  -- breaks optical purity without any D-LDH present
CLASSES: list[tuple[str, str, str, str]] = [
    ("L_LDH", "TARGET", "(ec:1.1.1.27)",
     "NAD-dependent L-lactate dehydrogenase; the producer enzyme"),
    ("D_LDH", "DISQUALIFY", "(ec:1.1.1.28)",
     "NAD-dependent D-lactate dehydrogenase; takes a batch off-spec"),

    # L-LDH fold decoy.
    ("MDH", "DECOY", "(ec:1.1.1.37)",
     "malate dehydrogenase; same fold and PFAM as L-LDH"),

    # D-LDH fold decoys, the 2-hydroxyacid dehydrogenase family.
    ("SERA", "DECOY", "(ec:1.1.1.95)",
     "D-3-phosphoglycerate dehydrogenase; serine biosynthesis"),
    ("GHR", "DECOY", "(ec:1.1.1.79) OR (ec:1.1.1.81) OR (ec:1.1.1.26) OR (ec:1.1.1.29)",
     "glyoxylate/hydroxypyruvate/glycerate reductases"),
    ("PDXB", "DECOY", "(ec:1.1.1.290)",
     "erythronate-4-phosphate dehydrogenase; PLP biosynthesis"),
    ("D2HA", "DECOY", "(ec:1.1.1.345) OR (ec:1.1.99.6) OR (ec:1.1.1.272)",
     "assorted D-2-hydroxyacid dehydrogenases"),
    ("FDH_NAD", "DECOY", "(ec:1.17.1.9)",
     "NAD-dependent formate dehydrogenase; 2-hydroxyacid DH family"),

    # Lactate consumers. 1.1.2.x / 1.1.5.x hand electrons to the respiratory chain,
    # so unlike 1.1.1.x they are effectively one-way and consume rather than produce.
    #
    # These are queried by gene name, not EC. UniProt has demoted lldD to EC 1.1.-.-
    # (checked: P33232, E. coli LldD), so `ec:1.1.2.3` now returns zero reviewed
    # bacterial entries and an EC-based consumer screen silently finds nothing.
    ("LLDD", "CONSUME", "(gene:lldD) OR (gene:lctD) OR (gene:lutA)",
     "L-lactate dehydrogenase, FMN/cytochrome-linked; eats L-lactate"),
    ("DLD", "CONSUME", "(gene:dld) OR (ec:1.1.5.12) OR (ec:1.1.2.4) OR (ec:1.1.2.5)",
     "D-lactate dehydrogenase, quinone/cytochrome-linked; eats D-lactate"),
    ("LOX", "CONSUME", "(ec:1.13.12.4) OR (gene:lctO) OR (gene:lox)",
     "lactate 2-monooxygenase / lactate oxidase"),

    # Optical purity without any D-LDH present.
    ("LARA", "PURITY", "(ec:5.1.2.1) OR (gene:larA)",
     "lactate racemase LarA; racemises pure L-lactate"),
    ("LAR_MAT", "PURITY",
     "(gene:larB) OR (gene:larC) OR (gene:larE) OR (gene:larC1) OR (gene:larC2)",
     "LarB/LarC/LarE nickel-pincer cofactor maturation; LarA is inactive without them"),
]

def fetch(query: str, fmt: str = "fasta", tries: int = 4) -> str:
    url = f"{STREAM}?query={query}&format={fmt}"
    for attempt in range(tries):
        proc = subprocess.run(
            ["curl", "-sS", "--max-time", "120", "--compressed", url],
            capture_output=True, text=True,
        )
        if proc.returncode == 0 and proc.stdout.strip():
            return proc.stdout
        time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"UniProt fetch failed: {query}\n{proc.stderr[:400]}")


def quote(s: str) -> str:
    from urllib.parse import quote as q
    return q(s, safe="")


def parse_fasta(text: str) -> list[tuple[str, str]]:
    records, header, seq = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:], []
        elif header is not None:
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    panel = OUT / "ldh_panel.faa"
    manifest = OUT / "panel_manifest.tsv"

    rows = []
    seen: dict[str, str] = {}          # accession -> class already assigned

    with panel.open("w") as fh:
        for label, role, ec_query, why in CLASSES:
            query = quote(f"({ec_query}) AND (reviewed:true) AND ({BACTERIA})")
            text = fetch(query)
            records = parse_fasta(text)

            kept = 0
            for header, seq in records:
                parts = header.split("|")
                acc = parts[1] if len(parts) > 2 else header.split()[0]

                # A protein carrying two of our ECs is ambiguous evidence, so it is
                # dropped rather than allowed to pull hits toward whichever class
                # happened to be fetched first.
                if acc in seen:
                    if seen[acc] != label:
                        seen[acc] = "AMBIGUOUS"
                    continue
                seen[acc] = label

                gene, org = "NA", "NA"
                for token in header.split():
                    if token.startswith("GN="):
                        gene = token[3:]
                if "OS=" in header:
                    org = header.split("OS=")[1].split("OX=")[0].strip().replace(" ", "_")

                fh.write(f">{label}|{acc}|{gene}|{org}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
                kept += 1

            rows.append((label, role, why))
            print(f"  {label:9s} {role:11s} {kept:4d}  {why}", file=sys.stderr)

    # Drop any sequence whose accession turned out ambiguous across classes.
    ambiguous = {a for a, v in seen.items() if v == "AMBIGUOUS"}
    if ambiguous:
        keep_lines, skip = [], False
        for line in panel.read_text().splitlines():
            if line.startswith(">"):
                skip = line[1:].split("|")[1] in ambiguous
            if not skip:
                keep_lines.append(line)
        panel.write_text("\n".join(keep_lines) + "\n")
        print(f"  dropped {len(ambiguous)} accessions with conflicting EC classes",
              file=sys.stderr)

    # Count after ambiguous cross-class accessions have been removed so the manifest
    # describes the FASTA that is actually searched.
    final_counts: dict[str, int] = {}
    for line in panel.read_text().splitlines():
        if line.startswith(">"):
            label = line[1:].split("|", 1)[0]
            final_counts[label] = final_counts.get(label, 0) + 1
    with manifest.open("w") as fh:
        fh.write("class\trole\tn_sequences\trationale\n")
        for label, role, why in rows:
            fh.write(f"{label}\t{role}\t{final_counts.get(label, 0)}\t{why}\n")

    n_final = sum(1 for line in panel.read_text().splitlines() if line.startswith(">"))
    print(f"\npanel: {n_final} sequences -> {panel}", file=sys.stderr)

    subprocess.run(
        ["diamond", "makedb", "--in", str(panel), "-d", str(OUT / "ldh_panel"),
         "--quiet"],
        check=True,
    )
    print(f"diamond db: {OUT / 'ldh_panel.dmnd'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
