# What the 301 simulations already say about supporters

This is not a new result. It is the existing `comets/out/` landscape read for the one
question the pool design turns on, which does not appear to have been extracted before.

**Design:** one producer (`P1_B3`) crossed against eight supporters, communities of size
1–5, two media (glycerol-minimal and glycerol + yeast extract), 301 simulations.

## No community beat the producer alone. Not one.

| Medium | `P1_B3` monoculture | multi-member communities | any that beat it |
|---|---:|---:|---:|
| glycerol-minimal | 0.2351 | 149 | **0** |
| glycerol + yeast extract | 0.2794 | 150 | **0** |

## The typical cost is not marginal, it is total

Median producer biomass, as a fraction of its own monoculture:

| Partners added | glycerol-minimal | glycerol + yeast extract |
|---:|---:|---:|
| 1 | 2.0% | 1.3% |
| 2 | 1.8% | 1.2% |
| 3 | 1.6% | 1.1% |
| 4 | 1.5% | 1.0% |

The median partner does not reduce the producer, it erases it. And **yeast extract makes
it worse**, because supplementation feeds the fast grower more than it feeds the
producer.

## The ranking is growth rate, and nothing else

Producer biomass in each pairing on glycerol-minimal, against each partner's µ on
glycerol from the model gate:

| Partner | µ glycerol | producer left | partner biomass |
|---|---:|---:|---:|
| `P3_C11` *Gluconobacter* | 0.127 | **42.4%** | 0.149 |
| `P3_B12` *Lactococcus* | 0.179 | 25.3% | 0.237 |
| `P3_G8` *Rhodococcus* | 0.217 | 16.0% | 0.396 |
| `P2_G9` *Leclercia*-like | ~0.55 | 2.0% | 1.046 |
| `P2_D5` *Pantoea* | ~0.55 | 2.0% | 1.107 |
| `P2_D12` *Rahnella* | ~0.55 | 1.9% | 1.100 |
| `P3_F8` *Ewingella* | ~0.55 | 1.8% | 1.109 |
| `P2_E7` *Chryseobacterium* | 0.604 | 1.7% | 1.203 |

The order is the growth-rate order. On one shared carbon source the faster grower takes
it, and the producer's share is set by the ratio of the two rates. Nothing in the
supporter's annotated provisioning capacity shows up at all.

This is the same mechanism `screen_d_guild.py` documented (producer biomass 0.257 → 0.0033
alongside four supporters). The landscape shows it is not specific to those four: it is
what happens with **every** partner tested, at every community size, on both media.

## The precursor-supply hypothesis is worth about 19%, and no strain can deliver it

The yeast-extract arm is the clean test, because it supplies amino acids directly instead
of via an organism, so it isolates the benefit from the delivery cost:

- `P1_B3` monoculture, glycerol-minimal: **0.2351**
- `P1_B3` monoculture, glycerol + yeast extract: **0.2794**

Relieving biosynthetic burden is worth **+18.8%** to the producer. That is a real effect
and it justifies the hypothesis being taken seriously.

But the cheapest partner-mediated version of that same benefit costs **57.6%** of producer
biomass (`P3_C11`, the best case), and the median costs 98%. A partner that could supply
amino acids in quantity has to grow to do it, and growing is how it takes the carbon.

**The conclusion is not "communities don't work". It is that this particular benefit is
one you buy, not one you co-culture for.** Yeast extract, or a defined amino-acid
supplement, delivers +19% with no competitor in the vessel.

## What this does and does not establish

It does **not** rule out community benefit generally. COMETS/FBA cannot represent the
mechanism most likely to produce it: overflow secretion is not growth-coupled, so at
maximum growth these models secrete exactly zero of the provisions they carry exchanges
for. A real supporter that leaks amino acids while growing is invisible to this
simulation. Crude-glycerol impurity detoxification and lag-time reduction — the LS30
precedent's actual advantages — are also outside what was simulated here, and both are
plausible reasons to keep partners in the design.

It does establish that **on the axis the pool was built to test — partner raises producer
yield on a shared carbon source — the answer is uniformly negative across 299 communities,
and the magnitude of the loss tracks growth-rate mismatch precisely.** Any supporter kept
in the pool should be kept for impurity tolerance or lag, tested on crude glycerol, and
should be growth-rate matched. Keeping a µ = 0.60 organism alongside a µ = 0.12 producer
has a predictable outcome and the prediction has now been made 299 times.
