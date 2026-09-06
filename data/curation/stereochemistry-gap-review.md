# Stereochemistry and reaction completeness

The full 18,177-equation selenium-forward model contains 119 detected substrate/product pairs with different, fully specified stereochemistry and 46 pairs involving unspecified stereochemistry. These counts classify **pairs within existing equations**, not independently characterized Cannabis enzymes. Full equations, all coproducts and original sources are retained.

Among the 3,004 target records lacking a producing equation, 118 have a producer-linked counterpart with the same connectivity, isotope composition, protonation and bond orders after removing only stereo information:

- 63 records have a fully specified stereoisomer counterpart to investigate.
- 51 require unspecified stereochemistry resolved before proposing conversion.
- Four have prior SDF/XML source conflicts; source identity takes priority.

Eighteen records have a fully specified counterpart already carrying a conditional CO₂ certificate. This is a prioritization signal, **not 18 new pathways**. The candidates include (-)-isoborneol, cis-THC, selected terpene stereoisomers and glycosides. Names must be checked against exact structures; agreement between XML and SDF does not independently prove that a name is correct. Some inventory entries are environmental/exposure compounds rather than established plant biosynthetic products.

## Choosing the missing chemistry

Do not insert a generic reversible stereoisomerization edge solely because it balances. Search for an exact racemase/epimerase/isomerase reaction or a source-backed multistep route with stereoselective formation of the target.

For example, the borneol/isoborneol lead may warrant checking oxidation to a shared carbonyl intermediate followed by stereoselective reduction. [EC 1.1.1.198](https://enzyme.expasy.org/EC/1.1.1.198) and [EC 1.1.1.227](https://enzyme.expasy.org/EC/1.1.1.227) specify particular borneol/camphor stereoisomers and nicotinamide cofactors. Those entries do **not** by themselves establish reduction to the exact (-)-isoborneol target. That specificity remains a source-search or experimental hypothesis, not a completed route.

For every proposed addition, retain exact substrate and product stereochemistry, all required cofactors and coproducts, direction evidence and species-level support. Unknown enzyme capability remains an enzyme gap. A different encoded structure or an unspecified stereocentre remains an identity issue, not a biochemical transformation by default.

The audit adds no equations, merges no identities and changes no coverage. Its complete source-linked queue is `data/reports/phase1-stereochemistry-queue.json`.
