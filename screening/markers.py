#!/usr/bin/env python3
"""Marker gene sets for the glycerol-to-lactate screen.

Separated from the screen so the sets can be reviewed and cited independently of
the scoring. Names follow eggNOG ``Preferred_name`` conventions and deliberately
carry aliases from both Gram-negative (E. coli) and Gram-positive (B. subtilis,
lactic acid bacteria) nomenclature, because a single-nomenclature set silently
misses whole clades — the v1 screen penalised ``poxB`` but not its Gram-positive
equivalent ``poxL``, so the Lactococcus candidate escaped a pyruvate-drain
penalty it should have taken.

Absence of a named gene is weak evidence of absent function. Some pathways have
no canonical gene at all (alanine is made by promiscuous transaminases), so those
are marked unreliable and must not drive an auxotrophy call on their own.
"""

from __future__ import annotations

# ---------------------------------------------------------------- carbon entry

# Glycerol kinase route: glycerol -> G3P -> DHAP.
GLYCEROL_G3P_ROUTE = {"glpK", "glpD", "glpA", "glpB", "glpC", "gldA_g3p"}
# Aquaglyceroporins. MIP-family channels pass glycerol *and* dihydroxyacetone and
# other small polyols, so glpF is evidence for DHA movement in either direction,
# not only glycerol entry. v1 recorded it and never scored it.
GLYCEROL_FACILITATOR = {"glpF", "glpV", "fps1"}

# Cytoplasmic NAD-dependent glycerol dehydrogenase: glycerol -> DHA *inside* the
# cell, which then needs an exporter to be useful to a partner.
GLYCEROL_DH_CYTOPLASMIC = {"gldA", "dhaD", "gldB"}
# Membrane-bound PQQ-dependent glycerol/polyol dehydrogenase (acetic acid
# bacteria). Active site faces the periplasm, so DHA appears extracellularly
# without any dedicated transporter. Mechanistically distinct from gldA and the
# distinction decides whether a DHA shuttle can work at all.
GLYCEROL_DH_PERIPLASMIC = {"sldA", "sldB", "gldM", "sldH"}
PQQ_BIOSYNTHESIS = {"pqqA", "pqqB", "pqqC", "pqqD", "pqqE", "pqqF"}

# DHA phosphorylation. Two unrelated families: the PTS-dependent three-subunit
# DhaKLM, and a standalone ATP-dependent DhaK common in lactic acid bacteria and
# Enterobacteriaceae. v1 required all three subunits, so every ATP-dependent
# receiver in the pool was excluded, and donors carrying a lone dhaK were
# *rewarded* for "not consuming DHA".
DHA_KINASE_PTS = {"dhaK", "dhaL", "dhaM"}
DHA_KINASE_ATP = {"dhaK", "dak", "dakA", "dhaK1"}

LOWER_GLYCOLYSIS = {
    "tpiA", "gapA", "gap", "pgk", "gpmA", "gpmI", "gpmB", "pgm", "eno",
    "pyk", "pykA", "pykF",
}

# ------------------------------------------------------------------- lactate

L_LDH_KOS = {"K00016"}
D_LDH_KOS = {"K03778"}
MALATE_DH_KO = "K00024"
# Quinone/cytochrome-linked lactate dehydrogenases oxidise lactate; they are
# consumption evidence, never production.
LACTATE_OXIDATION_KOS = {"K00101", "K00102"}
LACTATE_CONSUMPTION = {"lldD", "dld", "lctD", "lldE", "lldF", "lctO", "lox"}

# Lactate transport. v1 had no export marker at all, which is why two candidates
# whose genomes carry no lactate transporter were ranked as direct producers and
# only failed later at the GEM stage. LctP/LldP-family permeases are the main
# annotated carriers; they are catabolic in origin but are the best available
# evidence that lactate can cross the membrane.
LACTATE_TRANSPORT = {"lctP", "lldP", "larD", "lutP", "yvfH"}

# --------------------------------------------- aerobic pyruvate competition

# Under aeration this is the dominant competitor for pyruvate: pyruvate ->
# acetyl-CoA -> TCA -> CO2 rather than -> lactate. It sets the yield ceiling and
# v1 did not look for it.
PYRUVATE_DEHYDROGENASE = {"pdhA", "pdhB", "pdhC", "pdhD", "aceE", "aceF", "lpdA", "lpd"}
# Pyruvate oxidase: pyruvate + O2 -> acetyl-P + H2O2. Aerobic pyruvate drain and
# a peroxide source.
PYRUVATE_OXIDASE = {"poxB", "poxL", "spxB", "pox", "cidC"}
# NADH oxidase regenerates NAD+ directly with O2, competing with LDH for NADH.
# In Lactococcus this is what shifts aerated cultures off homolactic metabolism.
#
# Gene name alone is NOT sufficient here. The `noxE` hit in the Lactococcus candidate
# carried no KEGG KO, no EC number, and a description of a different enzyme
# ("3-hydroxyacyl-CoA dehydrogenase, NAD binding domain") from the broad COG0446
# family — the same class of name collision as `cat` matching "Cation efflux family".
# A call therefore requires KO or description corroboration; the species literature
# may still support the phenotype, but the annotation on its own does not.
# Description matching is NOT usable either: "NADH:flavin oxidoreductase / NADH
# oxidase family" is the Old Yellow Enzyme family name, and matching it picked up
# trimethylamine dehydrogenase (tmd) and N-ethylmaleimide reductase (namA, nemA) in
# five pool members. Only the specific water-forming KOs count. Across this pool that
# yields zero calls, so the aerobic drain rests on PDH and pyruvate oxidase, which are
# well supported. A real NADH oxidase phenotype must come from species literature.
NADH_OXIDASE = {"noxE", "nox", "noxA", "noxB"}
NADH_OXIDASE_KOS = {"K00359", "K17869"}

COMPETING_SINKS = {
    "adhE", "adhA", "adhB", "pta", "ackA", "ack", "pflB", "pflD",
    "frdA", "frdB", "frdC", "frdD", "alsS", "alsD", "budA", "budB", "budC",
    "ldhA_pdo",
}
PDO_SINK = {"dhaB", "dhaC", "dhaE", "dhaT", "pduC", "pduD", "pduE", "pduQ"}
METHYLGLYOXAL = {"mgsA", "gloA", "gloB", "gloC"}

# ------------------------------------------------- oxygen handling / peroxide

# "cat" is deliberately excluded: it matched "Cation efflux family" in the
# Lactococcus candidate and chloramphenicol acetyltransferase in the Bacillus one,
# turning a catalase-negative strain into a false positive.
CATALASE = {"katA", "katB", "katE", "katG", "katX", "katN"}
PEROXIDASE_OTHER = {"ahpC", "ahpF", "tpx", "npr", "bcp", "osmC"}
SUPEROXIDE_DISMUTASE = {"sodA", "sodB", "sodC"}
TERMINAL_OXIDASE = {
    "cydA", "cydB", "cydC", "cydD", "qoxA", "qoxB", "qoxC", "qoxD",
    "ctaC", "ctaD", "ctaE", "ctaF", "cyoA", "cyoB", "cyoC", "cyoD", "ccoN", "ccoO",
}

# ---------------------------------------------------------------- biosafety

# Names worth a manual look, NOT a virulence determination. Gene-name matching
# against eggNOG Preferred_name is unreliable for this purpose because virulence
# nomenclature collides heavily with unrelated housekeeping genes: `mecA` matched
# the Bacillus/Lactococcus ClpC protease adaptor (nothing to do with methicillin
# resistance), and `esxA`, `spa`, `plc`, `inlA` and `sak` hit broad protein
# families present in non-pathogens. Those are excluded. A real assessment needs a
# curated virulence-factor database (VFDB/abricate), which is not installed here —
# so the screen's job is to say "this needs review", not to clear or condemn.
VIRULENCE_REVIEW_NAMES = {
    "hlyA", "hlyB", "hlyC", "hlyD", "hla", "hlb", "hld",
    "lukS", "lukF", "lukD", "lukE",
    "blaZ", "vanA", "vanB", "tst",
}
# Genera where a species-level identity and a documented risk-group assignment
# must precede lab work.
RISK_GENERA = {
    "Staphylococcus", "Streptococcus", "Enterococcus", "Klebsiella", "Citrobacter",
    "Acinetobacter", "Pseudomonas", "Mycobacterium", "Clostridium", "Bacillus",
    "Listeria", "Serratia", "Vibrio", "Gemella", "Cedecea", "Kluyvera",
    # Enterobacteriaceae relatives of Klebsiella that include opportunistic
    # pathogens. Raoultella reached the producer shortlist unflagged without these.
    "Raoultella", "Enterobacter", "Cronobacter", "Providencia", "Proteus",
    "Morganella", "Hafnia", "Lelliottia",
}

# ------------------------------------------------- amino acid biosynthesis

# ``unreliable`` marks pathways whose absence cannot be read as auxotrophy
# because the reaction is carried by promiscuous or unnamed enzymes.
AMINO_ACID_PATHWAYS: dict[str, dict[str, object]] = {
    "histidine": {"genes": {f"his{x}" for x in "ABCDEFGHI"}},
    "tryptophan": {"genes": {f"trp{x}" for x in "ABCDEFG"}},
    "arginine": {"genes": {f"arg{x}" for x in "ABCDEFGHIJ"} | {"carA", "carB"}},
    # Split, because valine/isoleucine and leucine fail independently and a
    # medium has to name which one it supplies.
    "valine_isoleucine": {
        "genes": {"ilvB", "ilvC", "ilvD", "ilvE"},
        "alternatives": [{"ilvH", "ilvC", "ilvD", "ilvE"}, {"ilvN", "ilvC", "ilvD", "ilvE"}],
    },
    "leucine": {"genes": {"leuA", "leuB", "leuC", "leuD"}},
    "threonine": {"genes": {"thrB", "thrC", "hom"}, "alternatives": [{"thrA", "thrB", "thrC"}]},
    # Homoserine acylation is metA or metX; sulfur incorporation is metB+metC
    # (transsulfuration) or metY (direct sulfhydrylation). All four combinations occur.
    "methionine": {
        "genes": {"metA", "metB", "metC"},
        "alternatives": [{"metA", "metY"}, {"metX", "metB", "metC"}, {"metX", "metY"}],
    },
    # dapD/dapE (succinylase), dapH/dapL (acetylase) and ddh (dehydrogenase) are
    # three alternative mid-pathway routes, so only the shared committed steps are scored.
    "lysine": {"genes": {"dapA", "dapB", "dapF", "lysA"}},
    "chorismate_shikimate": {
        "genes": {"aroA", "aroB", "aroC", "aroD", "aroE", "aroF", "aroG", "aroH",
                  "aroK", "aroL", "aroQ"}
    },
    "phenylalanine": {"genes": {"pheA", "pheB"}, "unreliable": True},
    "tyrosine": {"genes": {"tyrA", "tyrB", "tyrC"}, "unreliable": True},
    "serine": {"genes": {"serA", "serB", "serC"}},
    "glycine": {"genes": {"glyA", "glyA1"}, "unreliable": True},
    # Split: needing cysteine and needing a reduced sulfur source are different
    # medium requirements, and lumping them made both unreadable.
    "cysteine": {"genes": {"cysE", "cysK"}, "alternatives": [{"cysE", "cysM"}]},
    "sulfate_assimilation": {
        "genes": {"cysD", "cysN", "cysH", "cysI"},
        "alternatives": [{"cysD", "cysN", "cysH", "cysJ"}],
    },
    "proline": {"genes": {"proA", "proB", "proC"}},
    # GS/GOGAT or GDH; a strain with either assimilates ammonium.
    "glutamate_glutamine": {
        "genes": {"glnA", "gltB", "gltD"},
        "alternatives": [{"glnA", "gdhA"}],
    },
    "asparagine": {"genes": {"asnA", "asnB", "asnS"}, "unreliable": True},
    "alanine": {"genes": {"alaA", "alaC", "avtA", "alr"}, "unreliable": True},
    "aspartate": {"genes": {"aspC", "aspB"}, "unreliable": True},
}

# ------------------------------------------- vitamin / cofactor biosynthesis

# Sets are deliberately small and restricted to *committed* steps, with mutually
# exclusive route variants declared as alternatives. Large union sets spanning
# alternative routes cannot reach the completeness threshold by construction, which
# silently converts prototrophs into auxotrophs: a 6-gene pyridoxal set spanning
# both routes called every strain a B6 auxotroph even though P1_B3/P1_B11 are
# complete on pdxS+pdxT and P2_D10/P3_C9 are complete on pdxA+pdxJ+pdxH. The same
# flaw hid B12 synthesis in the Gluconobacters, which carry 17 cob* genes.
VITAMIN_PATHWAYS: dict[str, dict[str, object]] = {
    # Pyrimidine + thiazole branches; thiH (E. coli) and thiO (Bacillus) are
    # alternative glycine/tyrosine-derived thiazole routes.
    "thiamine_B1": {
        "genes": {"thiC", "thiD", "thiE", "thiG", "thiH"},
        "alternatives": [{"thiC", "thiD", "thiE", "thiG", "thiO"},
                         {"thiC", "thiD", "thiE", "thiG", "thiS"}],
        "salvage": {"thiM", "thiK"},
    },
    "riboflavin_B2": {"genes": {"ribA", "ribB", "ribD", "ribE", "ribH"}},
    # Pantothenate synthesis proper. The CoA steps (coaA-E) are near-universal and
    # not an auxotrophy axis, so they no longer dilute this call.
    "pantothenate_B5": {"genes": {"panB", "panC", "panD"}, "salvage": {"panE"}},
    "pyridoxal_B6": {
        "genes": {"pdxA", "pdxJ"},                    # DXP-dependent
        "alternatives": [{"pdxS", "pdxT"}],           # DXP-independent
        "salvage": {"pdxK", "pdxY", "pdxH"},
    },
    "biotin_B7": {"genes": {"bioF", "bioA", "bioD", "bioB"}},
    "folate_B9": {"genes": {"folE", "folK", "folP", "folC"}},
    "niacin_NAD_B3": {"genes": {"nadA", "nadB", "nadC"}, "salvage": {"pncA", "pncB", "nadV"}},
    # Expensive vitamin and a classic paid-for cross-feed. Corrin ring synthesis
    # runs by an aerobic (cob*) or anaerobic (cbi*) route; they are alternatives,
    # never both.
    "cobalamin_B12": {
        "genes": {"cobG", "cobI", "cobJ", "cobM", "cobH", "cobL", "cobN"},
        "alternatives": [{"cbiD", "cbiE", "cbiT", "cbiC", "cbiA", "cbiK", "cbiG"}],
        "salvage": {"btuB", "btuF"},
    },
    # Respiratory quinone, satisfied by ANY of three unrelated routes: classical
    # menaquinone (men*), the futalosine pathway (mqn*), or ubiquinone (ubi*).
    # Scoring only men* called the Pseudomonas candidates menaquinone auxotrophs
    # when they carry 8-9 ubi* genes and simply use the other quinone — an
    # over-called requirement inflates the medium it is meant to minimise.
    "respiratory_quinone": {
        "genes": {"menA", "menB", "menC", "menD", "menE", "menF", "menH"},
        "alternatives": [
            {"mqnA", "mqnB", "mqnC", "mqnD", "mqnE"},
            {"ubiA", "ubiB", "ubiC", "ubiD", "ubiE", "ubiF", "ubiG", "ubiH", "ubiI", "ubiX"},
        ],
    },
    # PDH cofactor. lipA/lipB synthesise it; lplA only salvages it. A
    # salvage-only strain requires exogenous lipoate — which also means a
    # lipoate-free medium suppresses PDH and pushes pyruvate toward lactate.
    "lipoate": {"genes": {"lipA", "lipB", "lipM"}, "salvage": {"lplA", "lipL"}},
    # Heme gates catalase and the cytochrome chain. Lactic acid bacteria are
    # typically heme auxotrophs, which is why they are catalase-negative and can
    # only respire when heme is supplied — a documented community mechanism and a
    # concrete medium lever under aeration.
    "heme": {"genes": {"hemA", "hemB", "hemC", "hemD", "hemE", "hemH"}},
}

# Amino acid exporters. Thin and clade-biased by nature: most amino acid efflux
# carriers are uncharacterised, so a supporter's *export* capability cannot be
# established from annotation and must come from the GEM or the bench.
AMINO_ACID_EXPORTERS = {
    "lysE", "rhtA", "rhtB", "rhtC", "yddG", "alaE", "brnF", "brnE", "eamA",
    "argO", "yeaS", "leuE", "mtr", "tyrP",
}

NITROGEN_ASSIMILATION = {"glnA", "gltB", "gltD", "gdhA", "amtB", "nirB", "nasA"}

# Cobalamin is only a growth requirement if the genome carries a B12-dependent
# enzyme AND lacks its B12-independent counterpart. Absence of B12 *synthesis* is
# not a B12 *requirement*: every panel genome carries metE, the B12-independent
# methionine synthase, so calling the whole pool B12-auxotrophic would have added
# an expensive vitamin to the medium for no reason.
B12_DEPENDENT_ENZYMES = {
    "metH": {"metE"},            # methionine synthase, B12-dependent vs B12-independent
    "nrdJ": {"nrdA", "nrdB", "nrdD"},   # class II ribonucleotide reductase
    "nrdJ1": {"nrdA", "nrdB", "nrdD"},
    "nrdJ2": {"nrdA", "nrdB", "nrdD"},
    # Methylmalonyl-CoA mutase and diol dehydratase are B12-dependent but sit on
    # catabolic/anaplerotic routes that are rarely essential, so they are recorded
    # rather than treated as creating a requirement.
}
B12_DEPENDENT_NONESSENTIAL = {"mutA", "mutB", "scpA", "pduC", "pduD", "pduE"}


def all_marker_names() -> set[str]:
    """Every gene name the screen recognises, for the marker-hit report."""
    names: set[str] = set()
    for group in (
        GLYCEROL_G3P_ROUTE, GLYCEROL_FACILITATOR, GLYCEROL_DH_CYTOPLASMIC,
        GLYCEROL_DH_PERIPLASMIC, PQQ_BIOSYNTHESIS, DHA_KINASE_PTS, DHA_KINASE_ATP,
        LOWER_GLYCOLYSIS, LACTATE_CONSUMPTION, LACTATE_TRANSPORT,
        PYRUVATE_DEHYDROGENASE, PYRUVATE_OXIDASE, NADH_OXIDASE, COMPETING_SINKS,
        PDO_SINK, METHYLGLYOXAL, CATALASE, PEROXIDASE_OTHER, SUPEROXIDE_DISMUTASE,
        TERMINAL_OXIDASE, VIRULENCE_REVIEW_NAMES, AMINO_ACID_EXPORTERS, NITROGEN_ASSIMILATION,
    ):
        names |= group
    for definition in (*AMINO_ACID_PATHWAYS.values(), *VITAMIN_PATHWAYS.values()):
        names |= set(definition.get("genes", set()))
        names |= set(definition.get("salvage", set()))
    return names
