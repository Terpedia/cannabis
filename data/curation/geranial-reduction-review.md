# Geranial reduction: exact chemistry versus Cannabis biology

Primary source: [ACS Catalysis, DOI 10.1021/acscatal.1c05334](https://pubs.acs.org/doi/10.1021/acscatal.1c05334), PMCID PMC8787751. Results and Figure 3B distinguish geranial supplied by geraniol oxidation from commercial citral-mixture assays. GluER from *Gluconobacter oxydans* forms the S product in an NADPH-regenerated cascade. Reported cascade conversion is 95.3%; S enantiomeric excess is 99.2%. These are different measurements, neither establishes exclusive S production, and neither is a universal reaction yield.

The added elementary channel is exact (E)-geranial + NADPH + H+ → exact (S)-citronellal + NADP+. Charged cofactors are reused from the Terpedia reaction catalog, not inferred from names. All five participants already existed. Independent RDKit-based element/isotope/charge balance and E/S identity checks pass. The catalog cofactor reference is not evidence of geranial specificity. No Cannabis protein is assigned and no assay glucose is introduced as a carbon seed.

## Full-inventory model result

The alcohol-acetates parent has 2,738 covered inventory records / 2,735 exact structural certificates. The new forward-only channel produces two additional exact net certificates:

| Inventory record | Exact modeled target | Steps | Net CO₂ per target |
|---|---|---:|---:|
| CDB000081 | (S)-citronellol | 28 | 10 |
| CDB000585 | S-configured citronellyl acetate | 29 | 12 |

The result is 2,740 conditionally covered records / 2,737 distinct certificates out of the unchanged 6,220-record inventory. All parent certificates are replayed unchanged. There are 18,222 balanced equations, 3,226 balanced-participant records, 2,987 records without a net-producing equation, 458 solver-reported infeasible records and 35 external inventory records. The added reaction and both complete certificates retain every participant, coefficient and direction.

## Important medium and pathway limitations

Both new certificates still consume imported iron–sulfur clusters, Fe(II), hydrogen, ammonium, hydrogen phosphate and protons alongside CO₂. These are dependencies of the permissive model, not a validated nutrient recipe. CO₂ is its sole carbon-containing input, but that alone does not demonstrate photosynthetic or physiologically usable synthesis. Pre-existing recycled pools, permissive inherited directions, energy, compartments and transport remain unresolved. The separately running no-cluster/no-peroxide analysis uses an older selenium-forward parent; its results cannot be attributed to this new scenario.

Next discriminating tasks are to test the new routes under curated uptake restrictions, review inherited directions and energy requirements, and identify/assay Cannabis candidate ene-reductases against pure geranial with chiral product analysis and explicit NADPH/NADH controls. The bacterial assay supports a chemical channel; it does not establish plant biosynthesis.

## Follow-up: targeted uptake restriction

`phase1-recent-gains-medium.json` applies the no-cluster/no-peroxide net-uptake restriction to the full 18,222-equation geranial model, solving only the four gains since the C17 release. This is separate from the ongoing full-inventory selenium-parent restriction. All four have independently validated alternative certificates: CDB000110 (30 steps), CDB000695 (30), CDB000081 (26), and CDB000585 (27). Their respective net CO₂ requirements are 12, 6, 10 and 12 per target.

Thus the imported clusters in the previously selected witnesses are not indispensable to these targets in this model. This does not establish physiological cluster assembly: the alternative certificates still consume Fe(II), H₂ and sulfide, with ammonium or—for CDB000695—N₂ among other inputs. The N₂ dependency is not evidence of Cannabis nitrogen fixation. Uptake forms, energetics, physiological reaction direction and minimality remain unresolved. These four results must not be extrapolated to all 6,220 targets or added to permissive coverage as new gains.
