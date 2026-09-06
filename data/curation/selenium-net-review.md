# Selenium incorporation sensitivity — 2026-09-06

## Result

The complete 6,220-record inventory was tested with the three added equations:

1. Exact neutral KEGG R03601: O-acetyl-L-serine + hydrogen selenide → L-selenocysteine + acetic acid.
2. Explicit computed speciation: hydrogen selenide ⇌ hydrogenselenide + H⁺.
3. Explicit computed speciation: acetic acid ⇌ acetate + H⁺.

All three conserve elements, isotopes and charge. The two organic source endpoints match the exact model structures, including encoded stereochemistry. The neutral hydrogen selenide and acetic acid endpoints are added as distinct structures, not merged with their charged forms.

Source: https://www.kegg.jp/entry/R03601. Original reaction text and four MOL files are pinned with retrieval timestamps and SHA256 checksums in `data/raw/selenium-kegg/retrievals.json`. KEGG's NAME mentions hydrogen **sulfide**, whereas its DEFINITION and exact compound equation specify **selenide**. This conflict remains explicit; EC 2.5.1.47 is a source annotation, not a Cannabis enzyme assignment.

| Scenario | Conditional covered records | New records | No producing equation | Solver-infeasible | Exchange records |
| --- | ---: | ---: | ---: | ---: | ---: |
| Parent local-speciation model | 2,723 | — | 3,004 | 458 | 35 |
| New incorporation forward only | 2,724 | 1 | 3,004 | 457 | 35 |
| New incorporation reversible | 2,724 | 1 | 3,004 | 457 | 35 |

Both sensitivities retain all 2,720 prior exact structural certificates and add one certificate for **CDB004952, selenocysteine**. There are now 18,177 balanced equations. Original direction exclusions and the full external-species boundary are unchanged. The 35 exchange records in the table are historical inventory records, not the count of allowed external species. The reversed new incorporation reaction provides no further target coverage in this test.

## What the certificate does and does not show

The selected forward certificate has 11 directed steps and consumes exactly three CO₂ per exported selenocysteine. All internal organic compounds are regenerated or net produced; no organic carbon is imported. The explicit source reaction and both proton bridges appear in the certificate.

This is **not a physiological Cannabis pathway**. The inherited permissive model admits hydrogen peroxide, iron redox species, and carbon-free iron–sulfur clusters as exchanges. This certificate uses them. It also uses reverse orientations of inherited catalog chemistry, including cysteine/acetate to O-acetylserine and taurine carboxylation. Restricting the newly added reaction to its forward orientation does not make those inherited orientations physiological. Pre-existing acetate and other internal pools can be recycled; their startup synthesis is not proven. Energetics, compartments, transport, enzymes and in-vivo flux remain unresolved.

The exact step list, extents, complete equations, external consumption and net exports are in `phase1-selenium-forward-net.json`; the separate reversible report must not be silently substituted for it.

## Remaining gaps and next evidence

Selenomethionine and selenohomocysteine remain solver-infeasible in both new scenarios. The prior selenium-connectivity obstruction describes the **parent** model only; it must not be presented as a proof against the extended selenocysteine model. Recompute the necessary cut for the extended model before attributing the remaining failures to a specific missing entry reaction.

Prioritize exact-source reactions connecting selenocysteine to selenocystathionine/selenohomocysteine, while checking stereochemistry, coproducts and all-input supply. Keep sulfur-to-selenium substrate analogies as testable hypotheses unless direct source evidence exists. Cannabis-specific enzyme capability will require separate evidence and testing; neither an EC join nor this net certificate establishes it.
