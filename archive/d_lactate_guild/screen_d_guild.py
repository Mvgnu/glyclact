#!/usr/bin/env python3
"""Assemble a D-lactate guild from the full library: every member is a producer.

This replaces the producer/supporter split, which the community landscape falsified.
Non-producing supporters were selected for what they could provide -- amino acids,
vitamins, catalase -- and two results say that is the wrong axis:

* FBA never secretes the provisions. Every supporter carries all 20 amino-acid
  exchanges and several mmol/gDW/h of capacity, and at maximum growth secretes
  exactly zero. Export is not growth-coupled.
* Supporters competitively exclude the producer. Producer biomass falls from 0.257
  in monoculture to 0.0033 alongside four supporters, ~80-fold, because the
  Enterobacteriaceae run at mu ~0.55-0.60 on glycerol against 0.117-0.129 for the
  Bacillus producers. On one shared substrate the fast grower simply takes it.

A partner that eats the feedstock and returns nothing measurable is a cost. So the
design inverts to the LS30 precedent (PMC6387908): a crude-glycerol enrichment whose
membership was Enterobacter 57.97%, Escherichia/Shigella 39.25%, Klebsiella 0.67%,
Citrobacter 0.41% -- four glycerol fermenters, a *redundant* guild rather than a
division of labour. Its advantage was impurity tolerance and lag time, not yield.

Redundancy also fixes the kinetic problem for free. A guild drawn from one physiology
is growth-rate matched by construction, so no member is excluded 80-fold by another.

Isomer discipline, because getting this wrong already cost one round of analysis. The
second digit of the EC number is the electron acceptor, and it decides direction:

  L-LDH   EC 1.1.1.27  K00016          NAD, reversible      makes L-lactate
  D-LDH   EC 1.1.1.28  K03778          NAD, reversible      makes D-lactate
  lldD    EC 1.1.2.3   K00101          cytochrome           *eats* L-lactate
  dld     EC 1.1.5.12  K03777          quinone              *eats* D-lactate
          EC 1.1.2.4   K00102          cytochrome           *eats* D-lactate

An earlier pass filtered on "any lactate oxidase" and cut the bench from 61 to 4 by
rejecting D producers for carrying lldD, which cannot touch their product. 1.1.1.x is
NAD(P)-linked and runs either way; 1.1.2.x and 1.1.5.x hand electrons to the respiratory
chain and are effectively one-way.

dld is a tier here and not a gate, but not because it is harmless. It is quinone-linked, so
it needs a terminal acceptor to turn over and is thermodynamically idle in a strictly
fermentative regime -- that argument stands on its own. It is *not* true, as an earlier
version of this docstring claimed, that E. coli retains dld while producing D-lactate:
Mazumdar et al. (AEM 2010) deleted it in both strains (LA01 dpflB dfrdA ddld; LA02 dfrdA
dpta dadhE ddld), reaching 32 g/L from 40 g/L glycerol at 85% of the theoretical maximum
and >99.9% chiral purity. The precedent removes the gene rather than tolerating it.

It stays a tier only because gating on it leaves 1 of 45 strains, and because it is a
deletion target rather than a property of the isolate. Under aeration or with nitrate
present it will re-consume product.

Annotation traps this handles explicitly, each one found by censusing the library rather
than by assuming the vocabulary:

* **EC substring collisions.** "1.1.1.28" occurs inside "1.1.1.284" (adhC) and "1.1.1.27"
  inside "1.1.1.271" (rmd). EC fields are tokenised on commas and compared whole.
* **Gene names outliving their orthology.** eggNOG emits `gldA` carrying ko:K01990, an ABC
  transporter, and `dhaT` carrying ko:K19954. A name is therefore only believed when the
  protein has no KO and no EC at all -- if an ortholog was assigned and it is not ours,
  the name is a stale synonym rather than evidence.
* **L-LDH against malate dehydrogenase.** Both sit in the same fold with the same PFAM
  architecture (Ldh_1_N + Ldh_1_C), so PFAM cannot separate them; the disambiguator is
  EC 1.1.1.37 / K00024 on the same protein. D-LDH is a different fold entirely
  (2-Hacid_dh + 2-Hacid_dh_C), so a D call is structurally checkable.
* **Lactate racemase.** The plausible-looking KO range K22227-K22232 is not lar: K22229 is
  gluconate 2-dehydrogenase. The real marker is K22373 / EC 5.1.2.1, and it appears once
  in 384 genomes.
* **Glycerol against propanediol dehydratase.** EC 4.2.1.30 (dhaB, K06120) acts on
  glycerol; EC 4.2.1.28 (pduCDE, K01699) is annotated to 1,2-propanediol and only
  dehydrates glycerol as a side reaction. They are counted separately.

Writes results/d_guild/d_guild_screen.tsv.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_GENOME_ROOT = Path(
    "/Volumes/T7 Shield/STRAIN_LIBRARY_SEQ_PROCCSSED_DATA/assmebled_and_proccessed_genomes"
)
LOCAL_GENOME_ROOT = Path(__file__).parent / "candidate_genomes"

COLUMN = {"query": 0, "description": 7, "name": 8, "ec": 10, "ko": 11, "pfam": 20}


@dataclass(frozen=True)
class Marker:
    """One enzyme, with the evidence that confirms it and the evidence that refutes it.

    `veto_ec`/`veto_ko` carry the near neighbour that eggNOG confuses this call with,
    so a hit can be actively disqualified rather than merely unsupported.
    """

    label: str
    ko: frozenset[str] = frozenset()
    ec: frozenset[str] = frozenset()
    names: frozenset[str] = frozenset()
    pfam: frozenset[str] = frozenset()
    veto_ec: frozenset[str] = frozenset()
    veto_ko: frozenset[str] = frozenset()


def marker(label, ko=(), ec=(), names=(), pfam=(), veto_ec=(), veto_ko=()) -> Marker:
    return Marker(
        label,
        frozenset(ko),
        frozenset(ec),
        frozenset(names),
        frozenset(pfam),
        frozenset(veto_ec),
        frozenset(veto_ko),
    )


# The lactate node. Every gate in this screen resolves here.
LACTATE = [
    marker("d_ldh", ko={"K03778"}, ec={"1.1.1.28"}, names={"ldhA", "ldhD"},
           pfam={"2-Hacid_dh", "2-Hacid_dh_C"}),
    marker("l_ldh", ko={"K00016"}, ec={"1.1.1.27"}, names={"ldhL", "lctE"},
           pfam={"Ldh_1_N", "Ldh_1_C"}, veto_ec={"1.1.1.37"}, veto_ko={"K00024"}),
    marker("dld_eats_D", ko={"K00102", "K03777"}, ec={"1.1.2.4", "1.1.5.12"}, names={"dld"}),
    marker("lldD_eats_L", ko={"K00101"}, ec={"1.1.2.3"}, names={"lldD", "lctO"}),
    marker("lactate_transporter", ko={"K03303"}, names={"lctP", "lldP"}),
    # A racemase turns an optically pure broth into a racemic one, so it is a purity
    # gate even though it is all but absent from this library.
    marker("racemase", ko={"K22373"}, ec={"5.1.2.1"},
           names={"larA", "larB", "larC", "larD", "larE"}),
]

# Glycerol entry. Which route runs decides where the electrons land, and therefore
# whether lactate is redox-favourable at all: GlpD passes them to the quinone pool
# (respiratory, needs an acceptor), GldA passes them to NAD (fermentative).
GLYCEROL = [
    marker("glpF", ko={"K02440"}, names={"glpF"}),
    marker("glpK", ko={"K00864"}, ec={"2.7.1.30"}, names={"glpK"}),
    marker("glpD_quinone", ko={"K00111"}, ec={"1.1.5.3"}, names={"glpD", "glpA"}),
    marker("glpABC_anaerobic", ko={"K00112", "K00113"}, names={"glpB", "glpC"}),
    marker("gldA_nad", ko={"K00005"}, ec={"1.1.1.6"}, names={"gldA", "dhaD"}),
    marker("dhaK", ko={"K00863", "K05878", "K05879", "K05881"}, names={"dhaK", "dhaL", "dhaM"}),
]

# Products that compete for the same NADH lactate needs. LS30 made 27.77 g/L 1,3-PDO
# against 14.68 g/L lactate: on glycerol the reductive branch is the main rival, not
# acetate. Klebsiella, Citrobacter and Enterobacter are exactly the genera that carry it.
COMPETING_PRODUCTS = [
    marker("glycerol_dehydratase", ko={"K06120"}, ec={"4.2.1.30"}, names={"dhaB", "dhaB1"}),
    marker("propanediol_dehydratase", ko={"K01699", "K13919", "K13920"}, ec={"4.2.1.28"}),
    marker("pdo_oxidoreductase", ko={"K00086"}, ec={"1.1.1.202"}),
    marker("acetoin_synthase", ko={"K01652"}, ec={"2.2.1.6"}, names={"alsS", "budB"}),
    marker("acetoin_decarboxylase", ko={"K01575"}, ec={"4.1.1.5"}, names={"alsD", "budA"}),
    marker("butanediol_dh", ko={"K00004"}, ec={"1.1.1.4"}, names={"budC", "bdhA"}),
    marker("ethanol_adhE", ko={"K04072"}, names={"adhE"}),
    marker("acetate_pta", ko={"K00625", "K13788"}, ec={"2.3.1.8"}, names={"pta"}),
    marker("acetate_ackA", ko={"K00925"}, ec={"2.7.2.1"}, names={"ackA"}),
    marker("formate_pfl", ko={"K00656"}, ec={"2.3.1.54"}, names={"pflB"}),
    marker("succinate_frd", ko={"K00244", "K00245", "K00246", "K00247"},
           names={"frdA", "frdB", "frdC", "frdD"}),
]

# Crude glycerol is the feedstock: residual methanol, soap and free long-chain fatty
# acids, 50-200 g/L Na/K, and heavy metals. In a redundant guild this is where members
# stop being interchangeable -- they cover different impurities.
ROBUSTNESS = [
    marker("methanol_oxidation", ko={"K14028", "K14029", "K23995", "K17066"},
           names={"mxaF", "mxaI", "xoxF", "mdo"}),
    marker("formaldehyde_detox", ko={"K00121", "K01070", "K08093", "K08094", "K00148"},
           names={"frmA", "frmB", "adhC", "fghA", "hxlA", "hxlB"}),
    marker("lcfa_transport", ko={"K06076"}, names={"fadL"}),
    marker("lcfa_activation", ko={"K01897"}, ec={"6.2.1.3"}, names={"fadD"}),
    marker("beta_oxidation", ko={"K06445", "K00632", "K01782", "K01825"},
           names={"fadE", "fadA", "fadB"}),
    marker("osmoprotectant", ko={"K00108", "K00130", "K02168", "K02000", "K02001",
                                 "K02002", "K06718", "K00836", "K06720"},
           names={"betA", "betB", "betT", "proV", "proW", "proX", "ectA", "ectB", "ectC"}),
    marker("na_antiport", ko={"K03313", "K03314", "K05565"}, names={"nhaA", "nhaB", "mrpA"}),
    marker("metal_efflux", ko={"K15726", "K17686", "K01534", "K03325", "K00520"},
           names={"czcA", "copA", "zntA", "arsB", "merA"}),
    marker("catalase", ko={"K03781", "K03782"}, ec={"1.11.1.6"}, names={"katE", "katG"}),
]

ALL_MARKERS = LACTATE + GLYCEROL + COMPETING_PRODUCTS + ROBUSTNESS


def tokens(value: str, strip_prefix: str = "") -> set[str]:
    """Split an eggNOG multi-value field. Whole tokens only -- EC numbers nest."""
    out = set()
    for token in value.split(","):
        token = token.strip()
        if strip_prefix:
            token = token.removeprefix(strip_prefix)
        if token and token != "-":
            out.add(token)
    return out


def scan(annotation: Path) -> tuple[Counter[str], dict[str, list[str]]]:
    counts: Counter[str] = Counter()
    evidence: dict[str, list[str]] = {m.label: [] for m in ALL_MARKERS}

    with annotation.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= COLUMN["pfam"]:
                fields.extend(["-"] * (COLUMN["pfam"] + 1 - len(fields)))

            name = fields[COLUMN["name"]].strip()
            ecs = tokens(fields[COLUMN["ec"]])
            kos = tokens(fields[COLUMN["ko"]], "ko:")
            pfams = tokens(fields[COLUMN["pfam"]])

            for probe in ALL_MARKERS:
                if (probe.veto_ec & ecs) or (probe.veto_ko & kos):
                    counts[f"{probe.label}__vetoed"] += 1
                    continue

                by_ko = bool(probe.ko & kos)
                by_ec = bool(probe.ec & ecs)
                # An EC that is present and disagrees outranks the ortholog call:
                # eggNOG assigns KOs by best-hit, and EC comes from the curated mapping.
                if by_ko and ecs and probe.ec and not by_ec:
                    counts[f"{probe.label}__ec_conflict"] += 1
                    continue

                if by_ko or by_ec:
                    counts[probe.label] += 1
                    evidence[probe.label].append(name or fields[COLUMN["query"]])
                    if probe.pfam and probe.pfam <= pfams:
                        counts[f"{probe.label}__fold_confirmed"] += 1
                elif name and name in probe.names and not kos and not ecs:
                    # Names are believed only in the absence of any orthology call.
                    counts[probe.label] += 1
                    counts[f"{probe.label}__name_only"] += 1
                    evidence[probe.label].append(name)

    return counts, evidence


def provenance(name: str) -> str:
    """Not every genome in the library is a strain you can pull out of a freezer.

    81 of the 384 are GCF reference assemblies downloaded for comparison, and the MYb/JUb
    series is the C. elegans natural-microbiome collection. They are scanned and reported,
    because excluding them silently would hide what the screen actually ranked, but a pool
    can only be built from isolates the lab holds.
    """
    if name.startswith("GCF_"):
        return "reference_assembly"
    if name[:3] in {"P1_", "P2_", "P3_"}:
        return "lab_plate"
    if name.startswith(("MYb", "JUb", "BIGb", "CEent", "MSPm")):
        return "worm_collection"
    return "other"


def read_taxonomy(genome_dir: Path) -> tuple[str, float, str]:
    """MEGAN genus with an honest confidence.

    unique_taxa_percentage.txt normalises *within each rank row*, so a row holding one
    taxon always reads 100 -- which is how P1_B3 came to be labelled "Aneurinibacillus,
    100%" while its family row said Bacillaceae 99.98%. unique_taxa_stats.txt carries
    the real coverage, so that is what is reported, and a genus whose coverage falls
    well short of its family's is flagged rather than trusted.
    """
    stats = genome_dir / "diamond_megan" / "unique_taxa_stats.txt"
    if not stats.exists():
        return "", 0.0, "no_taxonomy"

    rows = [line.split() for line in stats.read_text(errors="replace").splitlines() if line.strip()]
    parsed = [(r[0], float(r[-1])) for r in rows if len(r) >= 3]
    if not parsed:
        return "", 0.0, "no_taxonomy"

    genus, genus_coverage = parsed[-1]
    family_coverage = parsed[-2][1] if len(parsed) >= 2 else genus_coverage
    if genus_coverage < 40.0:
        flag = "weak"
    elif family_coverage - genus_coverage > 25.0:
        flag = "genus_uncertain"
    else:
        flag = "clean"
    return genus, round(genus_coverage, 1), flag


def classify(counts: Counter[str]) -> dict[str, object]:
    """Gates, then tiers. Never one additive score.

    The last screen that summed weighted traits ranked the three strains already in the
    pool at the top, on generic stress genes. Breadth beat fitness for the role. Here a
    strain either belongs to the guild or does not, and everything else is a named tier
    so a trade-off stays visible instead of being averaged away.
    """
    d_ldh = counts["d_ldh"]
    l_ldh = counts["l_ldh"]
    glycerol_kinase = counts["glpK"] and counts["glpF"]
    glycerol_fermentative = counts["gldA_nad"] and counts["dhaK"]

    # Carrying D-LDH does not make a strain a producer. Lactate is an overflow product of a
    # fermentative pyruvate node, and an obligate respirer routes pyruvate through PDH to
    # CO2 instead. Without this gate the bench filled with Pseudomonas -- 14 of them, none
    # carrying pflB or adhE, all of which oxidise glycerol rather than spilling it. Every
    # Enterobacteriaceae genus here carries both, which is why the gate reproduces the LS30
    # membership rather than fighting it.
    fermentative = counts["formate_pfl"] > 0 or counts["ethanol_adhE"] > 0

    gates = {
        "makes_D_lactate": d_ldh > 0,
        "no_L_lactate_dh": l_ldh == 0,
        "no_racemase": counts["racemase"] == 0,
        "uses_glycerol": bool(glycerol_kinase or glycerol_fermentative),
        "ferments_pyruvate": fermentative,
    }
    failed = [name for name, passed in gates.items() if not passed]

    if glycerol_kinase and glycerol_fermentative:
        route = "both"
    elif glycerol_kinase:
        route = "glpK_respiratory"
    elif glycerol_fermentative:
        route = "gldA_fermentative"
    else:
        route = "none"

    pdo = counts["glycerol_dehydratase"] > 0 and counts["pdo_oxidoreductase"] > 0
    butanediol = counts["acetoin_synthase"] > 0 and counts["acetoin_decarboxylase"] > 0
    if pdo and butanediol:
        rival = "1,3-PDO + 2,3-BD"
    elif pdo:
        rival = "1,3-PDO"
    elif butanediol:
        rival = "2,3-BD"
    else:
        rival = "acids only"

    # Idle without a terminal acceptor, active under aeration. Reported so the choice
    # between a strict bench and a knockout stays the user's.
    reuptake = "aerobic risk" if counts["dld_eats_D"] else "none"

    impurities = [
        label for label, present in (
            ("methanol", counts["methanol_oxidation"] > 0 and counts["formaldehyde_detox"] > 0),
            ("soap_lcfa", counts["lcfa_activation"] > 0 and counts["beta_oxidation"] >= 2),
            ("salt", counts["osmoprotectant"] >= 3 and counts["na_antiport"] > 0),
            ("metals", counts["metal_efflux"] >= 2),
        ) if present
    ]

    return {
        "guild_member": int(not failed),
        "failed_gates": ",".join(failed) or "-",
        "d_ldh_copies": d_ldh,
        "d_ldh_fold_confirmed": counts["d_ldh__fold_confirmed"],
        "l_ldh_copies": l_ldh,
        "l_ldh_calls_rejected_as_mdh": counts["l_ldh__vetoed"],
        "reuptake_risk": reuptake,
        "eats_D_lactate": counts["dld_eats_D"],
        "eats_L_lactate_only": counts["lldD_eats_L"],
        "exports_lactate": counts["lactate_transporter"],
        "racemase": counts["racemase"],
        "glycerol_route": route,
        "ferments_pyruvate": int(fermentative),
        "competing_product": rival,
        "pdo_genes": counts["glycerol_dehydratase"] + counts["pdo_oxidoreductase"],
        "acetate_genes": counts["acetate_pta"] + counts["acetate_ackA"],
        "ethanol_genes": counts["ethanol_adhE"],
        "formate_lyase": counts["formate_pfl"],
        "succinate_genes": counts["succinate_frd"],
        "impurities_handled": ",".join(impurities) or "none",
        "impurity_breadth": len(impurities),
        "catalase": counts["catalase"],
    }


def main() -> None:
    root = Path(__file__).parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-root", type=Path, default=DEFAULT_GENOME_ROOT)
    parser.add_argument("--output", type=Path,
                        default=root / "results" / "d_guild" / "d_guild_screen.tsv")
    parser.add_argument("--top", type=int, default=25)
    options = parser.parse_args()

    genome_root = options.genome_root
    if not genome_root.is_dir():
        genome_root = LOCAL_GENOME_ROOT
        print(f"drive unavailable; scanning the local panel copies in {genome_root}")

    rows: list[dict[str, object]] = []
    directories = sorted(p for p in genome_root.iterdir() if p.is_dir())
    skipped = []
    for position, directory in enumerate(directories, start=1):
        annotation = directory / "eggnog.emapper.annotations"
        if not annotation.exists():
            skipped.append(directory.name)
            continue
        counts, _ = scan(annotation)
        genus, coverage, flag = read_taxonomy(directory)
        rows.append({
            "genome": directory.name,
            "source": provenance(directory.name),
            "genus": genus or "unassigned",
            "genus_coverage_percent": coverage,
            "taxonomy_flag": flag,
            **classify(counts),
        })
        if position % 50 == 0 or position == len(directories):
            print(f"scanned {position}/{len(directories)}", flush=True)

    if not rows:
        raise SystemExit(f"no annotations found under {genome_root}")

    options.output.parent.mkdir(parents=True, exist_ok=True)
    rows.sort(key=lambda r: (-int(r["guild_member"]),
                             -int(r["impurity_breadth"]),
                             r["competing_product"] != "acids only",
                             -int(r["exports_lactate"]),
                             str(r["genome"])))
    with options.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    guild = [r for r in rows if r["guild_member"]]
    bench = [r for r in guild if r["source"] != "reference_assembly"]
    print(f"\nscanned {len(rows)} genomes -> {options.output}")
    if skipped:
        print(f"{len(skipped)} directories carried no eggNOG annotation and were not scanned")

    print("\nattrition")
    stages = [
        ("all genomes scanned", rows),
        ("carry D-LDH", [r for r in rows if r["d_ldh_copies"]]),
        ("D-LDH and no L-LDH", [r for r in rows if r["d_ldh_copies"] and not r["l_ldh_copies"]]),
        ("also ferments pyruvate (pflB or adhE)",
         [r for r in rows if r["d_ldh_copies"] and not r["l_ldh_copies"]
          and r["ferments_pyruvate"]]),
        ("also glycerol-capable, no racemase (guild)", guild),
        ("guild, excluding reference assemblies", bench),
        ("bench that exports lactate", [r for r in bench if r["exports_lactate"]]),
        ("bench with no dld (strict, no knockout needed)",
         [r for r in bench if not r["eats_D_lactate"]]),
    ]
    for label, subset in stages:
        print(f"  {label:<48}{len(subset):>5}")

    print("\nguild by provenance")
    for source, n in Counter(str(r["source"]) for r in guild).most_common():
        print(f"  {source:<48}{n:>5}")

    print("\ncompeting reductive product across the bench")
    for rival, n in Counter(r["competing_product"] for r in bench).most_common():
        print(f"  {rival:<48}{n:>5}")

    print("\ngenus composition of the bench")
    for genus, n in Counter(str(r["genus"]) for r in bench).most_common(14):
        print(f"  {genus:<48}{n:>5}")

    print("\nbench ranked by impurity breadth, then clean product spectrum")
    print(f"{'genome':<10}{'genus':<17}{'route':<19}{'rival':<12}{'reuptake':<14}"
          f"{'impurities':<25}{'exp':>4}")
    for row in bench[: options.top]:
        print(f"{row['genome']:<10}{str(row['genus'])[:16]:<17}"
              f"{str(row['glycerol_route']):<19}{str(row['competing_product'])[:11]:<12}"
              f"{str(row['reuptake_risk']):<14}{str(row['impurities_handled'])[:24]:<25}"
              f"{row['exports_lactate']:>4}")


if __name__ == "__main__":
    main()
