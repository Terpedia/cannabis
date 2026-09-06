# Amino-phospholipid reaction gaps

This batch addresses exact PE, phosphatidylserine (PS), monomethyl-PE and dimethyl-PE structures in the 6,220-record inventory. It does not assign enzymes or relax encoded stereochemistry.

## Identity and speciation

The source-scaffold audit matches 162 records: 78 PE, 78 PS, three monomethyl-PE and three dimethyl-PE. Of these, 161 lacked producing equations in the published source-mapped-protonation model. One already had a conditional certificate.

The PE-family source forms are zwitterionic. Converting 84 neutral encoded structures to those source forms changes the location of protons and formal charges without changing the net formula or charge. These are explicitly labeled internal proton-relocation hypotheses. The 78 PS proposals instead include explicit H+ exchange. The two cases are not conflated. Neither is an exact ChEBI endpoint mapping, and the source reaction's identifier is evidence for its scaffold, not a claim that Rhea curates the proposed speciation equation.

Neutral and charged forms remain separate nodes. Isotopes, all heavy-atom bonds, and encoded stereochemistry are preserved. The audit rejects missing or inverted glycerol stereochemistry rather than filling it in. Its matching rule also excludes non-hydrocarbon or cyclic acyl tails. A non-match is a limitation of this bounded template set, not evidence that a compound cannot be synthesized.

## Synthesis proposals

| Source-forward template | Required inputs | Products | Instantiated equations |
| --- | --- | --- | ---: |
| [RHEA:32944](https://www.rhea-db.org/rhea/32944) | CDP-ethanolamine + exact diacylglycerol | PE + CMP + H+ | 78 |
| [RHEA:16914](https://www.rhea-db.org/rhea/16914) | Exact CDP-diacylglycerol + L-serine | PS + CMP + H+ | 78 |
| [RHEA:11165](https://www.rhea-db.org/rhea/11165) | Exact PE + S-adenosyl-L-methionine | Monomethyl-PE + S-adenosyl-L-homocysteine + H+ | 3 |
| [RHEA:32736](https://www.rhea-db.org/rhea/32736) | Exact monomethyl-PE + S-adenosyl-L-methionine | Dimethyl-PE + S-adenosyl-L-homocysteine + H+ | 3 |

These 162 synthesis equations accompany 162 separate speciation hypotheses. Every equation is checked for element, isotope and charge balance. Independent forward structural replay reconstructs every exact product from the generated reactants; atom-order reversal is tested separately. Synthesis is constrained to the source-written forward orientation as an explicit model assumption, not as a demonstrated physiological direction in Cannabis. Speciation is separately assumed reversible.

## Upstream supply

Fifty-six initial precursor structures were absent from the parent inventory: inputs to 51 PS targets and five methylated-PE targets. Backward expansion of all 165 distinct required input structures, including those already in the inventory, examines 278 structures and generates 177 precursor equations:

- 78 CDP-diacylglycerol synthesis proposals using CTP and exact phosphatidate ([RHEA:16230](https://www.rhea-db.org/rhea/16230));
- 78 sn-2 acylations and 15 sn-1 acylations with explicit acyl-CoA donors ([RHEA:19710](https://www.rhea-db.org/rhea/19710), [RHEA:15326](https://www.rhea-db.org/rhea/15326));
- three upstream PE syntheses and three first methylations.

All 56 missing initial precursors have proposed producing equations. The remaining 101 input structures occur in the parent inventory, but inventory presence is **not** treated as a supply certificate. The frontier remains explicitly recorded for the net-conversion model to evaluate. No carbon-containing donor, cofactor or lipid is added as an external seed.

## Interpretation and testing hypotheses

The completed amino-phospholipid sensitivity adds 374 distinct equations; 127 proposals join equations already represented. Its 17,789-equation model yields 96 independently replayed new exact net certificates: 66 PS, 28 PE, one monomethyl-PE and one dimethyl-PE. All previous 2,414 certified target records remain covered. Conditional coverage is now 2,510/6,220 (40.4%), while balanced participation covers 2,964 records.

Of the 161 formerly missing-producer records addressed here, 96 acquire certificates and 65 acquire producing equations but remain solver-reported infeasible. Thus 3,248 inventory records still lack producing equations, 427 are solver-reported infeasible, and 35 are exchange species. The existence of precursor structures and balanced reactions is not enough to close those remaining pathway gaps.

The speciation, synthesis and precursor reports alone add zero CO2 route claims. The separate full-inventory net report establishes conditional chemical feasibility under the unchanged CO2-only external carbon boundary. Such certificates permit regenerated pre-existing pools; they do not establish startup from zero pools, energy feasibility, compartment transport or Cannabis enzyme activity.

Discriminating future work includes resolving the unassigned glycerol configurations in excluded records, measuring relevant tissue/compartment pH and lipid speciation, testing exact-chain substrates rather than assuming generic lipid specificity, and comparing labeled serine, ethanolamine and methyl-donor incorporation against the competing synthesis routes. An assay or annotation must be linked to the exact tested reaction before it is counted as biological support.
