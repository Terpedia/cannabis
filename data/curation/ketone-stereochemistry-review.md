# Ketone-mediated stereochemistry sensitivity

This is a **reaction-class analogy**, not evidence of target-specific enzyme activity in Cannabis. Exact structures, rather than the labels below, define the targets. Source names may require further identity curation.

The generator changes a specified secondary alcohol to a ketone, retaining all other encoded stereochemistry, connectivity, charge and isotopes. A target and a specified stereoisomer counterpart must yield the same exact ketone. Unspecified stereochemistry and known source-identity conflicts are excluded. No direct, unrestricted stereoisomerization edge is added.

Each pair proposes explicit cofactor-balanced steps:

1. Partner alcohol + NAD⁻ → ketone + NADH²⁻ + H⁺.
2. Ketone + NADH²⁻ + H⁺ → target alcohol + NAD⁻.

The exact borneol/camphor equation in the Terpedia reaction catalog supplies reaction-class and cofactor bookkeeping evidence only. Its EC reference does **not** establish the substrate scope or stereoselectivity of an enzyme for another target. The report retains the complete source equation and directed source identifiers. Newly added equations are forward-only hypotheses; two exact equations already in the model retain their inherited directions.

## Full-input results

All 6,220 inventory records were assessed. Twenty balanced proposals produce 18 new equations and two exact existing-equation joins. All 2,721 previous structural certificates remain valid. Six additional exact structural certificates raise conditional record coverage from 2,724 to 2,730; this is **not confirmed biological coverage**.

| Record | Source label | Steps in saved certificate |
| --- | --- | ---: |
| CDB006165 | Maltose | 23 |
| CDB000071 | L-Arabinose | 12 |
| CDB000184 | (-)-Isoborneol | 28 |
| CDB000244 | L-Quebrachitol | 22 |
| CDB000467 | (S)-Ipsdienol | 28 |
| CDB004785 | Tetrahydrobiopterin | 34 |

Epicatechin (CDB004928) gains proposed producing chemistry but no full-input certificate. Missing-producer records decrease from 3,004 to 2,997, while solver-reported infeasible records increase from 457 to 458. Balanced participation is 3,215 records, distinct from pathway coverage. Thirty-five inventory records remain classified as external species rather than assessed synthesis targets.

Each saved positive certificate was replayed with exact rational stoichiometry, full element/isotope/charge balance, inherited direction bounds, nondepleting internal pools and CO₂ as the sole external carbon input. These witnesses retain the permissive 102-species exchange boundary. They are not results of the separate restricted-medium run and do not establish growth, zero-pool startup, energetics or physiological direction.

## Next discriminating work

- Resolve exact substrate stereochemistry against independently curated structures before trusting a source label.
- Seek substrate-specific reaction evidence for each oxidation and reduction separately. A known ketone is not evidence that the required stereoselective reduction exists.
- For the two certificates using only one newly added redox equation, inspect the existing ketone-producing route rather than claiming both proposed enzymes are required.
- Search Cannabis proteins for the relevant oxidoreductase capability, retaining homology, EC annotation and measured substrate activity as separate evidence classes. Assay the exact alcohol/ketone pair and distinguish stereoisomeric products; sequence similarity alone cannot settle stereoselectivity.
- Reassess these candidates under curated uptake and energy constraints. The current medium sensitivity uses the earlier selenium-forward model and must not be presented as validating this new scenario.

Reproducible artifacts: `phase1-ketone-stereo-hypotheses.json` and `phase1-ketone-stereo-net.json` in `data/reports/`, with lossless collection-row exports to Terpedia and matching `docs/data/` copies. Seven hypothesis/net and shared-builder regression tests passed, including replay of all previous and new certificates.
