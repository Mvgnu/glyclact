# Annotation-based candidate ranking (screen v2)

Priors for choosing models and experiments. They do not establish growth,
secretion, flux, yield or safety. Requirements are annotation-based and must
be confirmed against the GEMs; absence of a named gene is weak evidence.

## Direct glycerol-to-lactate producers (aerobic)

| Rank | Genome | Genus | % | Score | Isomer | Lactate direction | Export | Drain (PDH/POX/NOX) | Assembly | Requires | Roles |
|---:|---|---|---:|---:|---|---|---|---|---|---|---|
| 1 | `P3_A3` | Staphylococcus | 99.78 | 43.8 | L_and_D | producer_favoured | annotated | 4/1/0 | good | proline;serine;biotin_B7 | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 2 | `P1_B3` | Aneurinibacillus | 100.00 | 41.3 | L | producer_favoured | annotated | 4/0/0 | good | none | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 3 | `P1_B8` | Aneurinibacillus | 59.43 | 41.3 | L | producer_favoured | annotated | 4/0/0 | good | none | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 4 | `GCF_002979735.1_ASM297973v1_genomic` | unresolved | - | 40.8 | D | ambiguous | annotated | 2/1/0 | good | none | direct_producer,supporter,peroxide_scavenger |
| 5 | `GCF_002979955.1_ASM297995v1_genomic` | unresolved | - | 40.8 | D | ambiguous | annotated | 2/1/0 | good | none | direct_producer,supporter,peroxide_scavenger |
| 6 | `P1_B1` | Pseudomonas | 99.95 | 40.8 | D | ambiguous | annotated | 2/1/0 | good | none | direct_producer,supporter,peroxide_scavenger |
| 7 | `P1_B11` | Bacillus | 99.38 | 39.7 | L | producer_favoured | annotated | 4/0/0 | good | none | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 8 | `P1_F7` | Raoultella | 57.09 | 39.7 | L | producer_favoured | annotated | 4/0/0 | good | biotin_B7 | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 9 | `P3_F11` | Bacillus | 100.00 | 39.7 | L | producer_favoured | annotated | 4/0/0 | good | biotin_B7 | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 10 | `P2_D9` | Fictibacillus | 58.76 | 37.25 | L | producer_favoured | unannotated_gap | 4/0/0 | good | biotin_B7 | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 11 | `P3_A10` | Paenibacillus | 100.00 | 35.45 | L | producer_favoured | unannotated_gap | 4/1/0 | acceptable | none | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 12 | `P3_C5` | Paenibacillus | 100.00 | 35.3 | L | producer_favoured | unannotated_gap | 4/0/1 | good | biotin_B7 | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 13 | `GCF_002975115.1_ASM297511v1_genomic` | unresolved | - | 34.55 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 14 | `GCF_002975175.1_ASM297517v1_genomic` | unresolved | - | 34.55 | L | ambiguous | annotated | 4/0/0 | acceptable | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 15 | `GCF_002976395.1_ASM297639v1_genomic` | unresolved | - | 34.55 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 16 | `GCF_002976415.1_ASM297641v1_genomic` | unresolved | - | 34.55 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 17 | `GCF_002976445.1_ASM297644v1_genomic` | unresolved | - | 34.55 | L | ambiguous | annotated | 4/0/0 | acceptable | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 18 | `P1_A1` | Bacillus | 100.00 | 34.55 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 19 | `P1_A9` | Bacillaceae | 100.00 | 34.55 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 20 | `P1_B4` | Ureibacillus | 63.56 | 34.55 | L | ambiguous | annotated | 4/0/0 | acceptable | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 21 | `P1_C6` | Aneurinibacillus | 100.00 | 34.55 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 22 | `P1_H3` | Raoultella | 99.96 | 34.55 | L | ambiguous | annotated | 4/0/0 | acceptable | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 23 | `P2_F9` | Aneurinibacillus | 100.00 | 34.55 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 24 | `P3_B4` | Ureibacillus | 100.00 | 34.55 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 25 | `P3_B5` | Bacillaceae | 100.00 | 34.55 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |

## Glycerol-to-DHA donors

| Rank | Genome | Genus | % | Score | Isomer | Lactate direction | Export | Drain (PDH/POX/NOX) | Assembly | Requires | Roles |
|---:|---|---|---:|---:|---|---|---|---|---|---|---|
| 1 | `MYb596` | unresolved | - | 24 | none | consumer_favoured | annotated | 4/0/0 | no_quast | methionine;cobalamin_B12;niacin_NAD_B3(salvage_only)... | dha_donor,supporter |
| 2 | `MYb592` | unresolved | - | 21 | none | consumer_favoured | no_ldh | 4/0/0 | no_quast | methionine;respiratory_quinone | dha_donor,supporter,peroxide_scavenger |
| 3 | `MYb595` | unresolved | - | 21 | none | consumer_favoured | annotated | 4/0/0 | no_quast | methionine;pantothenate_B5;respiratory_quinone | dha_donor,supporter |
| 4 | `P3_C11` | Gluconobacter | 100.00 | 21 | none | consumer_favoured | no_ldh | 4/0/0 | good | methionine;pantothenate_B5;respiratory_quinone | dha_donor,supporter,peroxide_scavenger |
| 5 | `P3_C9` | Gluconobacter | 99.86 | 21 | none | consumer_favoured | no_ldh | 4/0/0 | good | methionine;pantothenate_B5;respiratory_quinone | dha_donor,supporter |
| 6 | `P3_H7` | Gluconobacter | 100.00 | 21 | none | consumer_favoured | no_ldh | 4/0/0 | good | methionine;pantothenate_B5;respiratory_quinone | dha_donor,supporter,peroxide_scavenger |
| 7 | `MYb375` | unresolved | - | 11 | D | ambiguous | annotated | 3/1/0 | no_quast | none | supporter,peroxide_scavenger |
| 8 | `MYb535` | unresolved | - | 11 | D | ambiguous | annotated | 3/1/0 | no_quast | none | supporter,peroxide_scavenger |
| 9 | `P1_G3` | Erwinia | 99.70 | 11 | D | consumer_favoured | annotated | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 10 | `CEent1` | unresolved | - | 10 | none | consumer_favoured | annotated | 3/1/0 | no_quast | pyridoxal_B6(salvage_only) | supporter,peroxide_scavenger |
| 11 | `P1_B10` | Rahnella | 57.26 | 10 | L_and_D | ambiguous | unannotated_gap | 3/1/0 | acceptable | none | supporter,peroxide_scavenger |
| 12 | `P1_C2` | Rouxiella | 50.90 | 10 | L_and_D | ambiguous | unannotated_gap | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 13 | `P1_C3` | Raoultella | 49.78 | 10 | L_and_D | ambiguous | unannotated_gap | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 14 | `P1_C4` | Rouxiella | 50.63 | 10 | L_and_D | ambiguous | unannotated_gap | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 15 | `P1_C5` | Rahnella | 52.11 | 10 | L_and_D | ambiguous | unannotated_gap | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 16 | `P2_E1` | Rahnella | 99.99 | 10 | L_and_D | ambiguous | unannotated_gap | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 17 | `P2_E2` | Rahnella | 99.99 | 10 | L_and_D | ambiguous | unannotated_gap | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 18 | `GCF_002979455.1_ASM297945v1_genomic` | unresolved | - | 8 | none | ambiguous | no_ldh | 5/1/0 | good | proline | supporter,peroxide_scavenger |
| 19 | `GCF_002979665.1_ASM297966v1_genomic` | unresolved | - | 8 | none | ambiguous | no_ldh | 5/1/0 | good | proline | supporter,peroxide_scavenger |
| 20 | `GCF_002979755.1_ASM297975v1_genomic` | unresolved | - | 8 | none | ambiguous | no_ldh | 5/1/0 | good | proline | supporter,peroxide_scavenger |
| 21 | `GCF_002980095.1_ASM298009v1_genomic` | unresolved | - | 8 | D | consumer_favoured | unannotated_gap | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 22 | `MYb121` | unresolved | - | 8 | D | consumer_favoured | unannotated_gap | 3/1/0 | no_quast | none | supporter,peroxide_scavenger |
| 23 | `MYb25` | unresolved | - | 8 | none | ambiguous | no_ldh | 5/1/0 | no_quast | proline | supporter,peroxide_scavenger |
| 24 | `MYb264` | unresolved | - | 8 | none | consumer_favoured | no_ldh | 5/1/0 | no_quast | proline | supporter,peroxide_scavenger |
| 25 | `MYb328` | unresolved | - | 8 | none | ambiguous | no_ldh | 5/1/0 | no_quast | proline | supporter,peroxide_scavenger |

## DHA-to-lactate receivers

| Rank | Genome | Genus | % | Score | Isomer | Lactate direction | Export | Drain (PDH/POX/NOX) | Assembly | Requires | Roles |
|---:|---|---|---:|---:|---|---|---|---|---|---|---|
| 1 | `P1_B3` | Aneurinibacillus | 100.00 | 43.4 | L | producer_favoured | annotated | 4/0/0 | good | none | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 2 | `P1_B8` | Aneurinibacillus | 59.43 | 43.4 | L | producer_favoured | annotated | 4/0/0 | good | none | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 3 | `P3_A3` | Staphylococcus | 99.78 | 42.5 | L_and_D | producer_favoured | annotated | 4/1/0 | good | proline;serine;biotin_B7 | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 4 | `P1_B11` | Bacillus | 99.38 | 41.8 | L | producer_favoured | annotated | 4/0/0 | good | none | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 5 | `P1_F7` | Raoultella | 57.09 | 41.8 | L | producer_favoured | annotated | 4/0/0 | good | biotin_B7 | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 6 | `P3_F11` | Bacillus | 100.00 | 41.8 | L | producer_favoured | annotated | 4/0/0 | good | biotin_B7 | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 7 | `P2_G9` | BIGb0383 | 99.88 | 37.5 | L_and_D | ambiguous | annotated | 3/0/0 | good | none | dha_receiver,supporter,peroxide_scavenger |
| 8 | `GCF_002975115.1_ASM297511v1_genomic` | unresolved | - | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 9 | `GCF_002975175.1_ASM297517v1_genomic` | unresolved | - | 37.0 | L | ambiguous | annotated | 4/0/0 | acceptable | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 10 | `GCF_002976395.1_ASM297639v1_genomic` | unresolved | - | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 11 | `GCF_002976415.1_ASM297641v1_genomic` | unresolved | - | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 12 | `GCF_002976445.1_ASM297644v1_genomic` | unresolved | - | 37.0 | L | ambiguous | annotated | 4/0/0 | acceptable | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 13 | `P1_A1` | Bacillus | 100.00 | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 14 | `P1_A9` | Bacillaceae | 100.00 | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 15 | `P1_B4` | Ureibacillus | 63.56 | 37.0 | L | ambiguous | annotated | 4/0/0 | acceptable | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 16 | `P1_C6` | Aneurinibacillus | 100.00 | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 17 | `P1_H3` | Raoultella | 99.96 | 37.0 | L | ambiguous | annotated | 4/0/0 | acceptable | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 18 | `P2_D9` | Fictibacillus | 58.76 | 37.0 | L | producer_favoured | unannotated_gap | 4/0/0 | good | biotin_B7 | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 19 | `P2_F9` | Aneurinibacillus | 100.00 | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 20 | `P3_B4` | Ureibacillus | 100.00 | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 21 | `P3_B5` | Bacillaceae | 100.00 | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 22 | `P3_C10` | Buttiauxella | 100.00 | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 23 | `P3_C4` | Paenibacillus | 99.77 | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 24 | `P3_C8` | Bacillaceae | 100.00 | 37.0 | L | ambiguous | annotated | 4/0/0 | good | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |
| 25 | `P3_F9` | Rahnella | 52.07 | 37.0 | L | ambiguous | annotated | 4/0/0 | fragmented | sulfate_assimilation | direct_producer,dha_receiver,supporter,peroxide_scavenger |

## Nutrient supporters

| Rank | Genome | Genus | % | Score | Isomer | Lactate direction | Export | Drain (PDH/POX/NOX) | Assembly | Requires | Roles |
|---:|---|---|---:|---:|---|---|---|---|---|---|---|
| 1 | `P2_C9` | Kluyvera | 61.98 | 66.0 | L_and_D | ambiguous | annotated | 3/1/0 | acceptable | none | supporter,peroxide_scavenger |
| 2 | `P3_A6` | Pseudomonas | 99.91 | 58.8 | D | consumer_favoured | annotated | 5/0/0 | good | none | supporter,peroxide_scavenger |
| 3 | `GCF_002980155.1_ASM298015v1_genomic` | unresolved | - | 57.8 | D | consumer_favoured | annotated | 2/0/0 | good | none | supporter,peroxide_scavenger |
| 4 | `MYb115` | unresolved | - | 57.8 | D | consumer_favoured | annotated | 2/0/0 | no_quast | none | supporter,peroxide_scavenger |
| 5 | `GCF_002979185.1_ASM297918v1_genomic` | unresolved | - | 56.8 | D | consumer_favoured | annotated | 2/0/0 | good | none | supporter,peroxide_scavenger |
| 6 | `GCF_002979635.1_ASM297963v1_genomic` | unresolved | - | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | dha_receiver,supporter,peroxide_scavenger |
| 7 | `GCF_002980055.1_ASM298005v1_genomic` | unresolved | - | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | dha_receiver,supporter,peroxide_scavenger |
| 8 | `GCF_002980135.1_ASM298013v1_genomic` | unresolved | - | 56.2 | D | consumer_favoured | unannotated_gap | 2/1/0 | good | none | supporter,peroxide_scavenger |
| 9 | `GCF_002980195.1_ASM298019v1_genomic` | unresolved | - | 56.2 | D | consumer_favoured | unannotated_gap | 2/1/0 | good | none | supporter,peroxide_scavenger |
| 10 | `GCF_009496595.1_ASM949659v1_genomic` | unresolved | - | 56.2 | D | consumer_favoured | annotated | 3/1/0 | acceptable | none | dha_receiver,supporter,peroxide_scavenger |
| 11 | `JUb66` | unresolved | - | 56.2 | L_and_D | ambiguous | annotated | 3/1/0 | no_quast | none | dha_receiver,supporter,peroxide_scavenger |
| 12 | `MYb186` | unresolved | - | 56.2 | L_and_D | ambiguous | annotated | 3/1/0 | no_quast | none | supporter,peroxide_scavenger |
| 13 | `P1_A10` | Raoultella | 96.87 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | dha_receiver,supporter,peroxide_scavenger |
| 14 | `P1_A12` | Citrobacter | 99.91 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 15 | `P1_A2` | Citrobacter | 99.90 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 16 | `P1_B12` | Citrobacter | 99.89 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 17 | `P1_B6` | Citrobacter | 99.93 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 18 | `P1_B7` | Citrobacter | 99.73 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 19 | `P1_C1` | Citrobacter | 99.95 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 20 | `P1_C10` | Citrobacter | 96.84 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 21 | `P1_C12` | Pantoea | 52.16 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 22 | `P1_D1` | Citrobacter | 99.89 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | good | none | supporter,peroxide_scavenger |
| 23 | `P1_D10` | Citrobacter | 95.92 | 56.2 | D | consumer_favoured | annotated | 3/1/0 | fragmented | none | supporter,peroxide_scavenger |
| 24 | `P1_D2` | Raoultella | 99.96 | 56.2 | L_and_D | ambiguous | annotated | 3/1/0 | good | none | dha_receiver,supporter,peroxide_scavenger |
| 25 | `P1_D4` | Raoultella | 99.96 | 56.2 | L_and_D | ambiguous | annotated | 3/1/0 | good | none | dha_receiver,supporter,peroxide_scavenger |

