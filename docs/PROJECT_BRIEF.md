# Project brief and decision rules

## Carbon accounting

Nutrient support alone does not establish glycerol-to-lactate conversion. The
carbon in lactate must originate from glycerol, either within a direct producer
or through a demonstrated extracellular intermediate such as DHA.

## Positive pathway evidence

### Glycerol-3-phosphate entry

`glpK` plus `glpD` or `glpABC`; `glpF` strengthens the call but is not mandatory
because alternative transport and passive entry can occur.

### DHA entry

`gldA`/`dhaD` plus `dhaKLM`. For a DHA shuttle, extracellular DHA accumulation
and growth of the receiving strain on DHA must be tested directly.

### Carbon flow to lactate

At least five of `tpiA`, `gapA`, `pgk`, `gpmA/gpmI`, `eno`, and `pyk`, plus an
NAD-dependent L- or D-lactate dehydrogenase. `lldD`, `dld`, and quinone-dependent
lactate dehydrogenases are lactate-consumption evidence, not production.

## Cautionary evidence

- `dhaBCE`/`pduCDE` and `dhaT`: glycerol to reuterin/1,3-propanediol.
- `mgsA-gloAB`: methylglyoxal route; interesting but stressful.
- `lldD`, `dld`, or `lctD`: possible lactate re-consumption.
- `adhE`, `pta-ackA`, `pflB`, `frdABCD`, `alsS-budABC`, or `poxB`:
  competing product capacity.

## Evidence gates

1. **Genome gate:** coherent pathways, plausible taxonomy, no obvious
   biosafety exclusion.
2. **Model gate:** mass-balanced model, growth on a biologically appropriate
   positive medium, no growth without a carbon source, and plausible glycerol
   uptake/product capacity.
3. **Reality gate:** measured single-strain substrate use and products agree
   directionally with the model. Models are corrected or rejected before
   community screening.
4. **Community gate:** simulated benefits require a named exchange mechanism
   and are tested in conditioned medium or co-culture.
5. **Carbon gate:** promising communities are confirmed using carbon balances
   and preferably labelled glycerol.

