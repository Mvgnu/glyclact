#!/usr/bin/env python3
"""Second-level signals: what distinguishes members of a deliberately redundant guild.

The upstream screens identify plausible L-lactate producers. This script asks what
functional diversity remains within a deliberately redundant, manually curated panel.

Redundancy was chosen on purpose -- it is what makes the guild growth-rate matched and what
LS30 actually demonstrated. But redundancy at the level of *lactate production* does not
mean redundancy everywhere, and the axes below are where members stop being interchangeable.

Three signals, none of them a ranking:

* **Gene-content distance.** Jaccard distance on KEGG ortholog profiles, which is a direct
  answer to "very different to the others". Reported as mean distance to the rest of the
  pool (an outlier score) and used to cluster the pool into representatives. This is the
  signal that exposes likely redundancy for manual review; it is not an automatic
  one-representative-per-cluster selector.
* **Biosynthetic self-sufficiency.** Amino acid and vitamin pathways, scored from KO
  presence. The point here is *not* provisioning -- FBA showed these strains never secrete
  amino acids at maximum growth, and non-producing provisioner partners were dropped for
  competitively excluding the producer. It is robustness: a prototroph runs in minimal
  medium on crude glycerol, an auxotroph needs the impurity profile to supply what it
  cannot make, which is exactly the variable that changes between feedstock batches.
  These KO-pathway calls are completeness proxies, not validated minimal-media phenotypes.
* **Accessory capability.** Cofactor-linked traits (B12, catalase, oxidase breadth) that
  determine whether a strain tolerates a different corner of the process.

Vitamin B12 is called from the *committed* corrinoid ring steps (cbiK/cobN and cobalamin
adenosyltransferase), not from the many cob-named genes that also serve siroheme and heme
biosynthesis -- a naive cob* grep marks almost every enterobacterium a B12 producer.

Fourteen producers is still a large design space, and the useful reduction is functional
rather than a second automated taxonomic gate.
ANI says how related two strains are; this says how differently they are equipped, which
is what decides whether a community of the two can do anything a monoculture cannot.

Writes results/pool_v5/guild_diversity.tsv and prints the functional clusters for review.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GENOMES_DIR = ROOT / "candidate_genomes"

# Committed steps, KO-first. A pathway is called complete when every listed step is found;
# partial counts are reported so a near-miss is visible rather than silently a zero.
AMINO_ACIDS: dict[str, list[set[str]]] = {
    # K11755 is the bifunctional hisIE that Enterobacteriaceae actually carry (44/45);
    # scoring only the split K01523/K01496 forms called every strain a His auxotroph.
    "His": [{"K00765"}, {"K01523", "K01496", "K11755"}, {"K00817"}, {"K00013"}],
    "Trp": [{"K01657", "K01658"}, {"K00766"}, {"K01695", "K01696", "K06001"}],
    "Phe": [{"K01850", "K04092", "K14170"}, {"K00832"}],
    # tyrA is K04517 (prephenate dehydrogenase) or K14187 (bifunctional with chorismate
    # mutase). The plausible-looking K00220/K24018/K15226 are absent from all 45.
    "Tyr": [{"K04517", "K14187"}, {"K00832"}],
    "Leu": [{"K01649"}, {"K01702", "K01703", "K01704"}, {"K00052"}],
    "Ile": [{"K01754"}, {"K01652", "K01653"}, {"K00053"}],
    "Val": [{"K01652", "K01653"}, {"K00053"}, {"K01687"}],
    "Lys": [{"K01714"}, {"K00674", "K00821"}, {"K01439", "K01438"}, {"K01586"}],
    "Thr": [{"K00928", "K12524"}, {"K00133"}, {"K01733"}],
    "Met": [{"K00651", "K00641"}, {"K01739", "K01740"}, {"K00548", "K00549"}],
    "Cys": [{"K00640"}, {"K01738"}],
    "Arg": [{"K00930", "K22478"}, {"K00145"}, {"K01755", "K14681"}],
    "Pro": [{"K00931"}, {"K00147"}, {"K00286"}],
    "Ser": [{"K00058"}, {"K00831"}, {"K01079"}],
    "Gly": [{"K00600"}],
    "Ala": [{"K00259", "K14260", "K19244"}],
    "Asp": [{"K00811", "K00812", "K00813"}],
    "Asn": [{"K01913", "K01953"}],
    "Glu": [{"K00261", "K00262", "K00265", "K00266"}],
    "Gln": [{"K01915"}],
}

VITAMINS: dict[str, list[set[str]]] = {
    "biotin_B7": [{"K00652"}, {"K00833"}, {"K01935", "K01012"}],
    "thiamine_B1": [{"K03147"}, {"K00788"}, {"K00941"}],
    "riboflavin_B2": [{"K01497"}, {"K00794"}, {"K00793"}],
    "pantothenate_B5": [{"K00606"}, {"K01918", "K13799"}],
    "pyridoxine_B6": [{"K03474", "K06215"}, {"K03473", "K00097"}],
    "folate_B9": [{"K01633"}, {"K00796"}, {"K11754", "K13939"}],
    "niacin_B3": [{"K03517"}, {"K00767"}],
    # De novo corrinoid ring construction only. The first attempt used K02302 and the
    # adenosyltransferase K00798 and called 44 of 45 B12 producers, which cannot be right:
    # E. coli makes no B12 de novo. K02302 is the cysG-family uroporphyrinogen-III
    # methyltransferase shared with siroheme, present in all 45, and K00798/btuR is a
    # *salvage* enzyme that adenosylates imported cobinamide. Neither implies de novo
    # synthesis. The committed steps are cobalt insertion and ring contraction (cbiK/cbiX,
    # cbiA) plus cobalamin-5'-phosphate synthase (cobS). The aerobic cobN route is absent
    # from the whole pool, so what these strains have is the anaerobic cbi pathway.
    "cobalamin_B12": [{"K02190", "K09882"}, {"K02224"}, {"K02233"}],
}

ACCESSORY: dict[str, set[str]] = {
    "catalase": {"K03781", "K03782"},
    "superoxide_dismutase": {"K04564", "K04565"},
    "cytochrome_bd_high_affinity": {"K00425", "K00426"},
    "nitrate_respiration": {"K00370", "K00371"},
    "fumarate_respiration": {"K00244", "K00245"},
    "hydrogenase": {"K06281", "K06282", "K18016", "K18017"},
    "b12_dependent_enzymes": {"K01847", "K00548"},
    "urease": {"K01428", "K01429", "K01430"},
    "cellulose_biofilm": {"K20542", "K20541"},
}


def assembly_qc(genome: str) -> dict[str, int]:
    """Assembly completeness, because it confounds every gene-content metric here.

    A fragmented or short assembly loses genes, which inflates gene-content distance and
    manufactures auxotrophies. Distinctness is only interpretable next to this: a strain
    that looks unique because it is missing a megabase is not unique, it is incomplete.
    """
    report = GENOMES_DIR / genome / "quast_output" / "report.tsv"
    if not report.exists():
        return {}
    values = {}
    for line in report.open(errors="replace"):
        if "\t" not in line:
            continue
        key, value = line.rstrip("\n").split("\t")[:2]
        try:
            values[key] = int(value)
        except ValueError:
            continue
    return {"contigs": values.get("# contigs", 0), "n50": values.get("N50", 0),
            "length": values.get("Total length", 0)}


def kegg_profile(annotation: Path) -> set[str]:
    kos: set[str] = set()
    with annotation.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            for token in fields[11].split(","):
                token = token.strip().removeprefix("ko:")
                if token and token != "-":
                    kos.add(token)
    return kos


def completeness(profile: set[str], pathway: list[set[str]]) -> tuple[int, int]:
    found = sum(1 for step in pathway if step & profile)
    return found, len(pathway)


def jaccard_distance(a: set[str], b: set[str]) -> float:
    union = a | b
    return 1.0 - len(a & b) / len(union) if union else 0.0


def functional_clusters(names: list[str], distance: dict[tuple[str, str], float],
                        threshold: float) -> list[list[str]]:
    """Cluster by single-linkage at a distance threshold, cheaply and deterministically.

    Nothing here needs scipy: the pool is 45 strains, and an explicit union-find keeps the
    result reproducible and easy to argue with.
    """
    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if distance[(a, b)] <= threshold:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
    groups: dict[str, list[str]] = {}
    for n in names:
        groups.setdefault(find(n), []).append(n)
    return sorted(groups.values(), key=lambda g: (-len(g), g[0]))


def main() -> None:
    root = ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path,
                        default=root / "candidate_genomes" / "manifest_pool_v5.tsv")
    parser.add_argument("--genomes", type=Path,
                        help="file with one genome id per line; overrides --manifest")
    parser.add_argument("--output", type=Path,
                        default=root / "results" / "pool_v5" / "guild_diversity.tsv")
    parser.add_argument("--cluster-threshold", type=float, default=0.10,
                        help="Jaccard distance below which two strains are near-duplicates.")
    options = parser.parse_args()

    if options.genomes:
        genomes = [ln.strip() for ln in options.genomes.read_text().splitlines() if ln.strip()]
    else:
        with options.manifest.open() as handle:
            genomes = [r["genome"] for r in csv.DictReader(handle, delimiter="\t")]
    profiles: dict[str, set[str]] = {}
    for genome in genomes:
        annotation = GENOMES_DIR / genome / "eggnog.emapper.annotations"
        if annotation.exists():
            profiles[genome] = kegg_profile(annotation)
    names = sorted(profiles)
    print(f"profiled {len(names)} genomes, {len(set().union(*profiles.values()))} distinct KOs")

    distance: dict[tuple[str, str], float] = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            d = jaccard_distance(profiles[a], profiles[b])
            distance[(a, b)] = distance[(b, a)] = d
    mean_distance = {
        n: sum(distance[(n, m)] for m in names if m != n) / max(1, len(names) - 1)
        for n in names
    }

    clusters = functional_clusters(names, distance, options.cluster_threshold)

    qc = {g: assembly_qc(g) for g in names}
    lengths = sorted(q["length"] for q in qc.values() if q.get("length"))
    median_length = lengths[len(lengths) // 2] if lengths else 0

    rows = []
    for genome in names:
        profile = profiles[genome]
        q = qc.get(genome, {})
        short = bool(median_length and q.get("length") and
                     q["length"] < 0.9 * median_length)
        fragmented = bool(q.get("contigs", 0) > 200)
        aa = {a: completeness(profile, p) for a, p in AMINO_ACIDS.items()}
        vit = {v: completeness(profile, p) for v, p in VITAMINS.items()}
        aa_complete = [a for a, (f, t) in aa.items() if f == t]
        aa_missing = [a for a, (f, t) in aa.items() if f < t]
        vit_complete = [v for v, (f, t) in vit.items() if f == t]
        accessory = [k for k, kos in ACCESSORY.items() if kos & profile]
        cluster = next(i for i, g in enumerate(clusters) if genome in g)
        rows.append({
            "genome": genome,
            "cluster": cluster,
            "cluster_size": len(clusters[cluster]),
            "distinctness": round(mean_distance[genome], 4),
            "distinctness_trustworthy": int(not (short or fragmented)),
            "assembly_contigs": q.get("contigs", 0),
            "assembly_length": q.get("length", 0),
            "assembly_caveat": ("short_assembly" if short else "") +
                               ("," if short and fragmented else "") +
                               ("fragmented" if fragmented else "") or "-",
            "kegg_orthologs": len(profile),
            "amino_acids_complete": len(aa_complete),
            "amino_acids_incomplete": ",".join(sorted(aa_missing)) or "-",
            "vitamins_complete": len(vit_complete),
            "vitamins_made": ",".join(sorted(vit_complete)) or "-",
            "makes_B12": int("cobalamin_B12" in vit_complete),
            "accessory_traits": ",".join(sorted(accessory)) or "-",
        })

    options.output.parent.mkdir(parents=True, exist_ok=True)
    rows.sort(key=lambda r: -float(r["distinctness"]))
    with options.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {options.output}\n")

    print(f"gene-content clustering at Jaccard <= {options.cluster_threshold}")
    print(f"  {len(names)} strains collapse to {len(clusters)} clusters")
    for i, group in enumerate(clusters):
        members = ", ".join(group)
        print(f"  cluster {i}: n={len(group):<3} {members}")

    print("\nmost distinct strains (mean gene-content distance to the rest)")
    print(f"  {'genome':<10}{'dist':>7}{'KOs':>7}{'AAs':>5}{'vitamins':>10}  accessory")
    for row in rows[:10]:
        print(f"  {row['genome']:<10}{row['distinctness']:>7.3f}"
              f"{row['kegg_orthologs']:>7}{row['amino_acids_complete']:>5}"
              f"{row['vitamins_complete']:>10}  {str(row['accessory_traits'])[:40]}")

    print("\nbiosynthetic self-sufficiency across the pool")
    full_aa = [r for r in rows if r["amino_acids_complete"] == len(AMINO_ACIDS)]
    print(f"  complete for all {len(AMINO_ACIDS)} amino acids : {len(full_aa)}/{len(rows)}")
    gaps = Counter()
    for r in rows:
        for a in str(r["amino_acids_incomplete"]).split(","):
            if a != "-":
                gaps[a] += 1
    if gaps:
        print("  most often incomplete: " + ", ".join(f"{a} ({n})" for a, n in gaps.most_common(6)))
    vitcount = Counter()
    for r in rows:
        for v in str(r["vitamins_made"]).split(","):
            if v != "-":
                vitcount[v] += 1
    print("\nvitamin pathways complete, by strain count")
    for v in VITAMINS:
        print(f"  {v:<20}{vitcount.get(v, 0):>4}/{len(rows)}")


if __name__ == "__main__":
    main()
