# Strain Identity and Lab Mapping

Bench reference for the 14-member fermentation-screen panel and the two purity-control
strains. Use the **assembly ID** in design matrices, Echo picklists, plate maps and analysis
files. Use the **lab marking** and **library well** to retrieve the physical strain.

Do not select a strain from biological name alone. Where a mapping is unresolved, verify it
against the original laboratory inventory before starting the experiment.

---

## 1. Fermentation-screen panel

| Assembly ID | Sequence-derived identity | Lab marking | Library well | Collection batch | Mapping status |
|---|---|---:|---|---|---|
| `P2_C2` | *Oerskovia merdavium* | `256` | `I8` | `Or` | CSV match |
| `P3_D2` | *Microbacterium* sp.; closest reference *M. forte* | `75` | `D1` | `Or` | CSV match |
| `P2_C10` | *Microbacterium* sp.; species unresolved | `269` | `J21` | `Or` | CSV match |
| `P3_G8` | *Rhodococcus erythropolis* | `294` | `I20` | `Or` | CSV match |
| `P3_C3` | *Lactococcus lactis* | `275` | `L1` | `Or` | CSV match |
| `P1_B11` | *Bacillus velezensis* | `41` | `B17` | `Or` | CSV match |
| `P1_F7` | *Bacillus altitudinis* group | `158` | `E16` | `Or` | CSV match |
| `P3_F11` | *Bacillus altitudinis* | `145` | `E21` | `Or` | CSV match |
| `P3_A9` | *Paenibacillus taichungensis* | `200` | `G6` | `Or` | CSV match |
| `P3_A8` | *Paenibacillus amylolyticus* | `196` | `G2` | `Or` | CSV match |
| `P3_C5` | *Paenibacillus odorifer* | `291` | `K13` | `Or` | CSV match |
| `P3_G9` | *Paenibacillus* sp.; closest reference *P. lupini* | `296` | `I22` | `Or` | CSV match |
| `P2_D9` | Bacillaceae sp.; closest reference *Niallia tiangongensis* | **unresolved** | **unresolved** | **unresolved** | absent from supplied CSV |
| `P2_B6` | *Pseudomonas* sp.; species unresolved | `264` | `I14` | `Or` | CSV match |

## 2. Purity-control strains

These two strains are not members of the 14-member acquisition panel. They are the `lar`
operon purity controls described in the fermentation-screen SOP and require the nickel
condition specified there.

| Assembly ID | Sequence-derived identity | Lab marking | Library well | Collection batch | Mapping status |
|---|---|---:|---|---|---|
| `P1_B3` | *Priestia aryabhattai* | `3` | `A3` | `Or` | CSV match |
| `P1_B8` | *Priestia megaterium* | **unresolved** | **unresolved** | **unresolved** | absent from supplied CSV |

---

## 3. Mapping provenance

Physical-stock fields were matched on `assembly_name` in:

```
/Users/magnusohle/Downloads/final_collection_excluding_dups_only_sequenced_with_assembly_name_FINAL.csv
```

Source columns:

| This document | CSV column |
|---|---|
| Assembly ID | `assembly_name` |
| Lab marking | `strain` (identical to `strain_Schullenberg_index_names` in all selected rows) |
| Library well | `Well.position` |
| Collection batch | `batch` |

Biological identities come from `panel/pool_v5_panel.yaml`; the CSV `Taxonomy` field is empty
for these records. Identity labels therefore describe the sequenced assemblies and do not
replace physical-stock verification.

### Unresolved entries

`P2_D9` and `P1_B8` do not occur in the supplied CSV. Their assembly folders and sequence data
exist in the project, but those files do not establish a physical lab marking or library well.
Leave both mappings unresolved until they are confirmed from the original inventory. Do not
infer them from a neighbouring assembly ID, species name or plate position.
