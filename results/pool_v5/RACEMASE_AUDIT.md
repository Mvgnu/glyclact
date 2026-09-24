# Lactate racemase audit

A racemase breaks optical purity with no D-LDH present anywhere, so D-LDH screening
cannot see it. `screen_d_guild.py` already looked for one — correctly rejecting the
plausible-looking KO range K22227–K22232 (K22229 is gluconate 2-dehydrogenase) in favour
of the real marker K22373 / EC 5.1.2.1 — and found it in 1 of 384 genomes.

That screen keyed on **LarA alone**. Two things were missed by doing so.

## LarA is inactive without its maturation operon — and the operon is the better probe

LarA needs a nickel-pincer cofactor (pyridinium-3,5-bisthiocarboxylic acid
mononucleotide, P2TMN) built by LarB, LarC and LarE. Screening the maturation genes as
well as LarA finds the machinery in **five** genomes rather than one, and LarA itself is
the divergent, hard-to-annotate member of the set — LarB/LarC/LarE are more conserved
and easier to detect. The operon is the more sensitive probe for the system.

## LarA is divergent enough to need an HMM, not a BLAST panel

Only three characterized bacterial LarA sequences exist in Swiss-Prot. The *Weissella*
LarA found here is 52% identical to the *L. plantarum* reference; anything more divergent
would fall below a homology search against three sequences. Searching the LarA
superfamily HMMs directly (PF09861 `Lar_N`, PF21113 `LarA_C`) is the sensitive test.

## Result

`results/pool_v5/lar_audit.tsv`

| Genome | Genus | larA | larB | larC | larE | Lar_N | LarA_C | Verdict |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `P2_A5` | *Lactiplantibacillus* | 99.8% | 99.6% | 100% | 99.3% | 252.7 | 237.3 | **full racemase** |
| `P3_G10` | *Weissella* | 52.1% | 51.9% | 41.5% | 59.2% | 232.8 | 200.5 | **full racemase** |
| `P1_B3` | (*Priestia*) | — | 55.4% | 47.7% | 38.8% | **39.5** | 0 | divergent LarA-superfamily protein inside the lar locus |
| `P1_B8` | (*Aneurinibacillus*) | — | 55.4% | 48.8% | 38.4% | **42.7** | 0 | same |
| `GCF_002979275.1` | — | — | 52.7% | 46.8% | 37.2% | 37.4 | 0 | same |

`P2_A5` (*L. plantarum*) and `P3_G10` (*Weissella*) carry complete, unambiguous lar
systems. Both also carry a D-LDH, so both were already excluded — the racemase does not
change their status, but it confirms the filter would have needed to catch them anyway.

## The finding that matters: `P1_B3`

**`P1_B3` is a producer in pool v4** — the keystone, per the archived
`archive/pool_v2_candidate_panel.yaml`. It
carries a contiguous, complete P2TMN cofactor operon:

```
BBFJKGGH_02803   larC1   nickel insertion protein
BBFJKGGH_02804   larC2   nickel insertion protein
BBFJKGGH_02805   larB    P2CMN synthase
BBFJKGGH_02806   larE    P2TMN synthase
BBFJKGGH_02808   ...     "hypothetical protein" -- Lar_N domain, 39.5 bits
```

All four maturation genes are reciprocal-best-hit confirmed at 38–55% identity. And the
"hypothetical protein" two genes along carries a **Lar_N domain** — the LarA superfamily
N-terminal domain. `P1_B8` has the identical arrangement.

### What this does and does not mean

It is **not** proof that `P1_B3` racemises lactate. Three things argue against:

- The Lar_N score is 39.5 bits against 252.7 for a real LarA, and the `LarA_C` domain is
  absent entirely.
- The locus is flanked by uronate metabolism genes (`uxuB`, uronate isomerase, an
  NagC-family repressor), not by lactate metabolism.
- The LarA superfamily is now known to be broad, with members acting as racemases and
  epimerases on several α-hydroxyacids other than lactate.

It is **not** dismissible either:

- A complete, dedicated nickel-pincer cofactor pathway is metabolically expensive and is
  not retained without a pincer-cofactor enzyme to serve. Something in `P1_B3` uses it.
- The candidate substrate-facing enzyme is annotated "hypothetical protein", which means
  no annotation-based screen of any kind would have flagged this strain.
- The genus label on this strain has already been wrong once — the archived panel YAML
  records that its "*Aneurinibacillus*, 100% genus support" call was a MEGAN
  normalisation artifact and ANI puts it in *Priestia aryabhattai*.

### Recommendation

This is an empirical question, and it is cheap to answer. `P1_B3` is already going into
culture; **run chiral HPLC on its monoculture supernatant** rather than reasoning further
from sequence. One injection settles it.

Until it is settled, `P1_B3` should not be the sole producer in any arm whose readout is
optical purity. It is exactly the hole that was warned about — a strain that can take a
batch off-spec with no D-LDH anywhere in its genome — and it was found in the one
producer the pool was built around.
