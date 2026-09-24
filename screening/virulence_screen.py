#!/usr/bin/env python3
"""Screen guild candidates for virulence and resistance by homology, not by genus.

"Citrobacter is an opportunistic pathogen" is a statement about a genus, and the genus is
not the thing that gets cultured. 48 of 61 candidates sat in genera with clinical
relatives, which is enough to stall a project and not enough to justify dropping a strain.
This asks the genome instead.

There is a direct analogue to ANI here, and it is the standard one: search the predicted
proteome against curated reference sets with DIAMOND and call a gene present above an
identity/coverage threshold. That is what ABRicate does. Two references are used:

* **VFDB core (setA)** -- 4,759 *experimentally verified* virulence factor proteins. The
  full setB adds predicted ones, which would bury real signal in homology noise.
* **NCBI AMRFinderPlus AMRProt** -- 9,998 resistance, stress and virulence proteins.

The two answer different questions and neither replaces ANI. ANI says what the organism is,
which is what maps onto a risk-group listing (TRBA 466 in Germany); this says what it
carries. A Pantoea agglomerans with no toxin and no acquired beta-lactamase is a different
proposition from a Klebsiella with a carbapenemase, even though both are "Enterobacteriaceae".

Reading the output requires manual biological judgment. The automated tier is a review
queue, not a biosafety verdict or a pool-selection gate. Every environmental Gram-negative
carries flagella, adhesins, siderophores and biofilm genes, and VFDB catalogues all of them
because they are virulence factors *in a pathogen's context*. Counting total VFDB hits would
rank the most ordinary soil bacterium as the most dangerous. So hits are split by VFDB class
and only three influence the review-priority heuristic:

    Exotoxin (VFC0235)                  a secreted toxin is a toxin anywhere
    Invasion (VFC0083)                  entering non-phagocytic cells is not environmental
    Effector delivery system (VFC0086)  reported, but see below

Effector delivery is the largest class in VFDB and the least specific: it holds both the
T3SS injectisome, which is a genuine pathogenicity marker, and the T6SS, which is an
antibacterial weapon carried by a great many harmless environmental Pseudomonas. It is
reported with its VF names attached rather than folded into the tier, because the name is
what distinguishes them.

Resistance is tiered on what actually constrains lab and disposal practice: carbapenemases,
ESBLs and colistin resistance. Arsenic, mercury, copper and quaternary-ammonium entries are
environmental adaptations and are counted separately from clinical classes.

Every hit prints the genes behind it. This project has already produced false positives
from name matching alone -- `cat` matching a cation efflux pump, `mecA` matching a ClpC
adaptor, `hlyC` and `vanB` raising flags that never survived inspection -- so a tier is
never reported without the evidence that produced it. Even so, complete-locus and
high-identity heuristics still make context-dependent or incorrect calls, so the final
inclusion/exclusion decision remains manual.

Writes results/screening/virulence_screen.tsv and results/screening/virulence_hits.tsv.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GENOME_ROOT = ROOT / "candidate_genomes"

# VFDB header: >VFG037176(gb|WP_001081735) (plc1) phospholipase C
#              [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii ACICU]
VFDB_HEADER = re.compile(
    r"^(?P<id>\S+)\s+\((?P<gene>[^)]*)\)\s+(?P<product>.*?)\s+"
    r"\[(?P<vf>.*?)\s+\((?P<vf_id>VF\d+)\)\s+-\s+(?P<vf_class>.*?)\s+\((?P<class_id>VFC\d+)\)\]"
    r"(?:\s+\[(?P<organism>.*?)\])?\s*$"
)

# A single hit does not make a pathogen. The first version of this screen treated the whole
# "Invasion" class as decisive and flagged 43 of 85 strains -- driven by ompA, the major
# outer membrane protein of every enterobacterium, which VFDB lists because it mediates
# blood-brain-barrier crossing in E. coli K1. ibeB, ibeC and aslA behave the same way.
#
# The fix is locus-level evidence: a multi-gene virulence system is only called when
# several of its genes are found together, so a structural protein borrowed by one
# pathogenic lineage cannot raise a flag on its own. Exotoxins are exempt because a
# secreted toxin is a single-gene determinant and is decisive by itself.
HIGH_CONSEQUENCE_SYSTEMS = (
    "K1 capsule", "K2 capsule", "TTSS", "T3SS", "SPI-1", "SPI-2", "Mxi-Spa", "Ysc",
    "LEE encoded TTSS", "Inv/Spa", "Psc", "Yersiniabactin", "Aerobactin", "RmpA",
)
LOCUS_GENE_MINIMUM = 3
# These loci are lineage-specific, so a decisive call demands near-identity rather than
# ordinary orthology. Three strains hit "K1 capsule" on gene count alone: one E. coli with
# all seven kps genes at 96-99%, and two Pantoea carrying only the kpsD/E/M export
# apparatus at 83-87%. The export machinery is shared by many capsule types; the identity
# floor is what separates the real locus from a conserved polysaccharide exporter.
HIGH_CONSEQUENCE_MIN_IDENTITY = 90.0

# T6SS is the largest signal in this bench and the least alarming: it is an antibacterial
# weapon carried by most environmental Pseudomonas and enterobacteria, not a determinant of
# human disease. It is reported by name and never tiered.
ANTIBACTERIAL_SYSTEMS = ("T6SS", "HSI-", "SCI ")

# AMRFinder's AMRProt doubles as the point-mutation reference set, so a plain homology
# search hits gyrA, parC, rpsL, folP and murA in every genome. Those carry no class, which
# is what keeps them out below. The `_rel` suffix is AMRFinder's own marker for proteins
# related to a resistance determinant but *not* conferring resistance -- an unfiltered run
# called qseB_rel/qseC_rel a colistin-resistance strain.
NON_RESISTANCE_SUFFIX = "_rel"

# Acquired determinants that would genuinely change how a strain is handled.
CRITICAL_FAMILIES = (
    "blaKPC", "blaNDM", "blaVIM", "blaIMP", "blaOXA-48", "blaOXA-23", "blaCTX-M",
    "blaGES", "blaSPM", "blaSIM", "mcr-", "vanA", "vanB", "cfiA",
)
# Chromosomal cephalosporinases native to the very genera on this bench: blaCMY is the
# Citrobacter freundii AmpC, blaACT the Enterobacter one, blaEC the E. coli one. Finding
# them is a statement about species identity, not about acquired resistance.
INTRINSIC_FAMILIES = (
    "blaCMY", "blaACT", "blaEC", "blaBUT", "blaYOC", "blaTER", "blaCFE", "blaDHA",
    "blaMIR", "blaFOX", "blaMOX", "blaPDC", "blaOKP", "blaLEN", "blaSHV-1", "ampC",
    "qnrB", "qnrE", "fosA", "oqxA", "oqxB", "acrB", "mdtM", "kdeA", "emrD",
)
ENVIRONMENTAL_AMR = {"ARSENIC", "MERCURY", "COPPER", "QUATERNARY_AMMONIUM", "EFFLUX",
                     "NICKEL", "TELLURIUM", "SILVER", "ACID", "BIOCIDE",
                     "COPPER/SILVER", "GOLD", "ZINC", "CADMIUM", "CHROMATE", "CHROMIUM",
                     "COPPER/NICKEL", "NICKEL/COBALT", "COBALT", "LEAD", "IRON"}
# AMRFinder files these under virulence rather than resistance.
AMR_VIRULENCE = {"STX1", "STX2", "INTIMIN", "HEAT_LABILE_ENTEROTOXIN", "HEAT_STABLE_ENTEROTOXIN"}


def run(command: list[str]) -> None:
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        sys.exit(f"command failed: {' '.join(command)}\n{process.stderr[-2000:]}")


def build_database(fasta: Path, database: Path, threads: int) -> None:
    if database.with_suffix(".dmnd").exists():
        return
    print(f"building {database.name}.dmnd")
    run(["diamond", "makedb", "--in", str(fasta), "-d", str(database),
         "--threads", str(threads), "--quiet"])


def vfdb_index(fasta: Path) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    with fasta.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            match = VFDB_HEADER.match(line[1:].rstrip("\n"))
            if not match:
                continue
            index[match.group("id")] = {
                "gene": match.group("gene"),
                "vf": match.group("vf"),
                "vf_class": match.group("vf_class"),
                "organism": match.group("organism") or "",
            }
    return index


def amr_index(fasta: Path) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    with fasta.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\n")
            fields = header.split("|")
            if len(fields) < 10:
                continue
            index[header.split()[0]] = {
                "gene": fields[3],
                "family": fields[4],
                "amr_class": fields[8],
                "product": fields[9].replace("_", " "),
            }
    return index


def selected_genomes(screen: Path | None, genome_list: Path | None) -> list[str]:
    """The strains to screen, from a plain id list or from a screen table.

    The id list is the current path: biosafety is a property of a strain, not of the
    screen that happened to select it, so this no longer requires the D-lactate guild
    table it was originally written against (now under archive/).
    """
    if genome_list is not None:
        chosen = [ln.strip() for ln in genome_list.read_text().splitlines() if ln.strip()]
        if not chosen:
            sys.exit(f"no genome ids in {genome_list}")
        return chosen

    with screen.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    chosen = [r["genome"] for r in rows]
    if not chosen:
        sys.exit(f"no genomes selected from {screen}")
    return chosen


def combined_proteome(genomes: list[str], genome_root: Path, target: Path) -> int:
    """One FASTA for all candidates, headers stamped with the genome they came from.

    DIAMOND is started once rather than per genome: startup dominates the runtime at this
    database size, and a single pass keeps memory flat.
    """
    written, missing = 0, []
    with target.open("w") as out:
        for genome in genomes:
            # Drive layout first, then the staged local copy that prepare_candidates.py
            # writes. The local copy is what lets this run with the library unmounted,
            # which is the same fallback identify_strains.py needs.
            faa = sorted((genome_root / genome / "prokka").glob("*.faa"))
            if not faa:
                faa = sorted((genome_root / genome).glob(f"{genome}.faa"))
            if not faa:
                missing.append(genome)
                continue
            for line in faa[0].read_text(errors="replace").splitlines():
                if line.startswith(">"):
                    out.write(f">{genome}@{line[1:].split()[0]}\n")
                else:
                    out.write(line + "\n")
            written += 1
    if missing:
        print(f"no prokka proteins for {len(missing)} genomes: {', '.join(missing[:6])}"
              + (" ..." if len(missing) > 6 else ""))
    return written


def search(query: Path, database: Path, output: Path, threads: int, identity: float) -> None:
    print(f"searching {database.name}")
    run(["diamond", "blastp", "-q", str(query), "-d", str(database), "-o", str(output),
         "-f", "6", "qseqid", "sseqid", "pident", "length", "slen", "evalue", "bitscore",
         "--id", str(identity), "--max-target-seqs", "1", "--evalue", "1e-10",
         # -b/-c bound the block size: the default sizing assumes far more RAM than this
         # machine has and would swap.
         "-b", "0.4", "-c", "1", "--threads", str(threads), "--quiet"])


def parse_hits(path: Path, index: dict[str, dict[str, str]], min_identity: float,
               min_coverage: float, kind: str) -> list[dict[str, object]]:
    hits = []
    seen: set[tuple[str, str]] = set()
    with path.open() as handle:
        for line in handle:
            query, subject, pident, length, slen, evalue, bitscore = line.rstrip("\n").split("\t")
            coverage = int(length) / int(slen)
            if float(pident) < min_identity or coverage < min_coverage:
                continue
            meta = index.get(subject)
            if not meta:
                continue
            genome, protein = query.split("@", 1)
            key = (genome, meta["gene"])
            if key in seen:
                continue
            seen.add(key)
            hits.append({
                "genome": genome,
                "protein": protein,
                "database": kind,
                "gene": meta["gene"],
                "identity_percent": round(float(pident), 1),
                "reference_coverage": round(coverage, 2),
                **{k: v for k, v in meta.items() if k != "gene"},
            })
    return hits


def classify(genome: str, hits: list[dict[str, object]]) -> dict[str, object]:
    """Tiers from decisive evidence, with the genes that produced them attached."""
    vf = [h for h in hits if h["database"] == "VFDB"]
    amr = [h for h in hits if h["database"] == "AMR"
           and not str(h["gene"]).endswith(NON_RESISTANCE_SUFFIX)]

    by_class = Counter(str(h.get("vf_class", "")) for h in vf)
    genes_per_system: dict[str, set[str]] = defaultdict(set)
    confident_genes_per_system: dict[str, set[str]] = defaultdict(set)
    for hit in vf:
        genes_per_system[str(hit["vf"])].add(str(hit["gene"]))
        if float(hit["identity_percent"]) >= HIGH_CONSEQUENCE_MIN_IDENTITY:
            confident_genes_per_system[str(hit["vf"])].add(str(hit["gene"]))

    high_consequence = sorted(
        system for system, genes in confident_genes_per_system.items()
        if len(genes) >= LOCUS_GENE_MINIMUM
        and any(marker.lower() in system.lower() for marker in HIGH_CONSEQUENCE_SYSTEMS)
    )
    antibacterial = sorted(
        system for system in genes_per_system
        if any(marker.lower() in system.lower() for marker in ANTIBACTERIAL_SYSTEMS)
    )
    exotoxins = [h for h in vf if h.get("vf_class") == "Exotoxin"]

    classed = [h for h in amr if h.get("amr_class")]
    environmental = [h for h in classed if h.get("amr_class") in ENVIRONMENTAL_AMR]
    amr_virulence = [h for h in classed if h.get("amr_class") in AMR_VIRULENCE]
    resistance = [h for h in classed if h not in environmental and h not in amr_virulence]

    critical = [h for h in resistance
                if any(str(h["gene"]).startswith(f) or str(h.get("family", "")).startswith(f)
                       for f in CRITICAL_FAMILIES)]
    intrinsic = [h for h in resistance
                 if any(str(h["gene"]).startswith(f) or str(h.get("family", "")).startswith(f)
                        for f in INTRINSIC_FAMILIES)]
    acquired = [h for h in resistance if h not in critical and h not in intrinsic]

    if exotoxins or amr_virulence or critical or high_consequence:
        tier = "review_required"
    elif acquired:
        tier = "elevated"
    else:
        tier = "background"

    reasons = []
    if exotoxins:
        reasons.append(f"exotoxin:{'/'.join(sorted({str(h['gene']) for h in exotoxins})[:4])}")
    if high_consequence:
        reasons.append(f"virulence_locus:{'/'.join(high_consequence[:3])}")
    if amr_virulence:
        reasons.append(
            f"toxin:{'/'.join(sorted({str(h['gene']) for h in amr_virulence})[:3])}")
    if critical:
        reasons.append(
            f"critical_resistance:{'/'.join(sorted({str(h['gene']) for h in critical})[:3])}")
    if acquired:
        reasons.append(f"acquired:{'/'.join(sorted({str(h['gene']) for h in acquired})[:4])}")

    return {
        "genome": genome,
        "biosafety_tier": tier,
        "tier_reasons": "; ".join(reasons) or "-",
        "exotoxin_genes": ",".join(sorted({str(h["gene"]) for h in exotoxins})) or "-",
        "high_consequence_loci": ",".join(high_consequence) or "-",
        "acquired_resistance_genes": ",".join(sorted({str(h["gene"]) for h in acquired})) or "-",
        "intrinsic_resistance_genes": ",".join(sorted({str(h["gene"]) for h in intrinsic})) or "-",
        "antibacterial_secretion": ",".join(antibacterial) or "-",
        "environmental_resistance_hits": len(environmental),
        "adherence_hits": by_class.get("Adherence", 0),
        "motility_hits": by_class.get("Motility", 0),
        "siderophore_hits": by_class.get("Nutritional/Metabolic factor", 0),
        "immune_modulation_hits": by_class.get("Immune modulation", 0),
        "invasion_class_hits": by_class.get("Invasion", 0),
        "total_vfdb_hits": len(vf),
    }


def main() -> None:
    root = ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-root", type=Path, default=DEFAULT_GENOME_ROOT)
    parser.add_argument("--screen", type=Path,
                        default=root / "results" / "pool_v5" / "screening_candidates.tsv",
                        help="table with a genome column; every row is screened")
    parser.add_argument("--genomes", type=Path,
                        help="file with one genome id per line; overrides --screen")
    parser.add_argument("--references", type=Path, default=root / "references")
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "screening")
    parser.add_argument("--work-dir", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-identity", type=float, default=80.0)
    parser.add_argument("--min-coverage", type=float, default=0.8)
    options = parser.parse_args()

    if not options.genome_root.is_dir():
        sys.exit(f"genome root not available: {options.genome_root}")

    vfdb_fasta = options.references / "vfdb_setA.fas"
    amr_fasta = options.references / "AMRProt.fa"
    for path in (vfdb_fasta, amr_fasta):
        if not path.exists():
            sys.exit(f"missing reference database: {path}")

    work = options.work_dir or (options.references / "diamond_work")
    work.mkdir(parents=True, exist_ok=True)

    genomes = selected_genomes(options.screen, options.genomes)
    print(f"{len(genomes)} candidate genomes selected from {options.screen.name}")

    query = work / "candidate_proteins.faa"
    scanned = combined_proteome(genomes, options.genome_root, query)
    print(f"combined proteomes of {scanned} genomes")

    build_database(vfdb_fasta, options.references / "vfdb_setA", options.threads)
    build_database(amr_fasta, options.references / "AMRProt", options.threads)

    vf_out, amr_out = work / "vfdb_hits.tsv", work / "amr_hits.tsv"
    search(query, options.references / "vfdb_setA", vf_out, options.threads, options.min_identity)
    search(query, options.references / "AMRProt", amr_out, options.threads, options.min_identity)

    hits = parse_hits(vf_out, vfdb_index(vfdb_fasta), options.min_identity,
                      options.min_coverage, "VFDB")
    hits += parse_hits(amr_out, amr_index(amr_fasta), options.min_identity,
                       options.min_coverage, "AMR")

    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for hit in hits:
        grouped[str(hit["genome"])].append(hit)

    rows = [classify(genome, grouped.get(genome, [])) for genome in sorted(genomes)]
    tier_rank = {"review_required": 0, "elevated": 1, "background": 2}
    rows.sort(key=lambda r: (tier_rank[str(r["biosafety_tier"])], -int(r["total_vfdb_hits"])))

    options.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = options.output_dir / "virulence_screen.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    hits_path = options.output_dir / "virulence_hits.tsv"
    columns = ["genome", "database", "gene", "identity_percent", "reference_coverage",
               "vf", "vf_class", "organism", "family", "amr_class", "product", "protein"]
    with hits_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore", restval="")
        writer.writeheader()
        writer.writerows(sorted(hits, key=lambda h: (str(h["genome"]), str(h["database"]),
                                                     str(h["gene"]))))

    print(f"\n{len(hits)} calls at >={options.min_identity:.0f}% identity over "
          f">={options.min_coverage:.0%} of the reference protein")
    print(f"wrote {summary_path}")
    print(f"wrote {hits_path}")

    print("\nbiosafety tiers")
    for tier, n in Counter(str(r["biosafety_tier"]) for r in rows).most_common():
        print(f"  {tier:<28}{n:>5}")

    flagged = [r for r in rows if r["biosafety_tier"] != "background"]
    if flagged:
        print("\nflagged strains, with the evidence")
        print(f"{'genome':<10}{'tier':<18}{'why':<60}")
        for row in flagged:
            print(f"{row['genome']:<10}{str(row['biosafety_tier']):<18}"
                  f"{str(row['tier_reasons'])[:59]:<60}")
    else:
        print("\nthe automated review heuristic found no high-priority determinant; "
              "manual hit review is still required")

    print("\nwhat the bench does carry, and why none of it tiers")
    totals = Counter()
    for row in rows:
        for system in str(row["antibacterial_secretion"]).split(","):
            if system != "-":
                totals[system] += 1
    for system, n in totals.most_common(6):
        print(f"  {system + ' (antibacterial T6SS)':<44}{n:>4} strains")
    intrinsic = Counter()
    for row in rows:
        for gene in str(row["intrinsic_resistance_genes"]).split(","):
            if gene != "-":
                intrinsic[gene.split("-")[0]] += 1
    for gene, n in intrinsic.most_common(6):
        print(f"  {gene + ' (species-typical, chromosomal)':<44}{n:>4} strains")
    for label, key in (("adherence", "adherence_hits"), ("motility", "motility_hits"),
                       ("siderophore", "siderophore_hits"),
                       ("ompA/ibe-type invasion homologs", "invasion_class_hits")):
        values = sorted(int(r[key]) for r in rows)
        print(f"  {label + ' (ubiquitous)':<44}median {values[len(values) // 2]:>3},"
              f" max {values[-1]:>3}")


if __name__ == "__main__":
    main()
