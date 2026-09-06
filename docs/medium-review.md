# Non-carbon inputs and the defined-medium workstream

[Current directed map](net.html?scenario=selenium&target=CDB004952) · [Exact medium audit](data/medium-inventory.json)

The current selenium-forward scenario has **2,724 conditional covered records out of 6,220**, represented by 2,721 exact structural certificates. The new selenocysteine certificate uses the incorporation reaction only forward. Allowing its reverse gives no additional coverage. This is a chemistry sensitivity, not demonstrated Cannabis metabolism.

## What is currently supplied?

The model permits **102 external species: CO₂ and 101 carbon-free species**. The saved certificates consume **CO₂ plus 22 distinct carbon-free inputs**. The 35 inventory entries classified as exchange species are target records—not the number of allowed medium components.

These are inputs used by particular selected routes, not proven essential nutrients:

| Net input | Exact structures consuming it |
| --- | ---: |
| Hydrogen phosphate | 629 |
| Sulfite | 17 |
| Nitrite | 1 |
| Hydrogen diphosphate | 1,260 |
| Hydrogen phosphite | 42 |
| Proton | 2,652 |
| Reduced 4Fe–4S reactive cluster | 2,503 |
| Ammonium | 413 |
| Reduced 2Fe–2S reactive cluster | 2,545 |
| Sulfate | 2 |
| Triphosphate | 11 |
| Chloride | 1 |
| Iron(II) | 2,523 |
| Hydrogen peroxide | 2,614 |
| Dinitrogen | 67 |
| Hydrogensulfide | 21 |
| Hydrogenselenide | 1 |
| Water | 4 |
| Oxidized 2Fe–2S reactive cluster | 32 |
| Phosphorous acid | 1 |
| Dihydrogen | 134 |
| Copper(I) | 29 |

Counts overlap because a route consumes multiple inputs. Exact charged structures, normalized amounts and all dependent target IDs are in the audit. Display names do not merge protonation states.

## Why this is not yet a minimum defined medium

Imported iron–sulfur clusters and peroxide expose a permissive redox boundary. A separate full-inventory test is evaluating blocked uptake of those species while allowing their production, regeneration and disposal. No completed result from that running test is claimed here.

Some model entries are reactive parts of proteins, not free nutrient molecules. [Photosystem I](https://www.rhea-db.org/rhea/30407) requires light and protein-bound plastocyanin/ferredoxin. [Photosystem II](https://www.rhea-db.org/rhea/36359) also requires photons. [Reviewed light requirements](data/light-reaction-requirements.json) are retained separately and attached to an edge when that reaction is present. Neither reaction appears in the current selected certificate set; these annotations are **not yet enforced energy constraints**.

Zero net mineral consumption does not establish dispensability: current certificates omit biomass dilution and many catalytic requirements. [Plant-nutrition guidance](https://extension.umd.edu/resource/garden-fertilizer-basics) identifies mineral requirements independently of this model. Dinitrogen use here is not evidence of direct nitrogen assimilation by Cannabis. Selenium demand for selenium-containing targets is separate from ordinary plant nutrient essentiality.

The next steps are to curate admissible uptake forms, minimize a shared input set, and replay complete pathways under that boundary. A feasible input set, an inclusion-minimal set and a proven minimum-cardinality set are different results. Light, starting pools, growth, compartments, transport, pH, concentrations and counterions remain outstanding. This is not a growth-medium recipe.

## Selenocysteine evidence

[Source structure audit](data/selenium-source-identity.json) · [Forward-only result](data/selenium-forward-net.json) · [Reversible comparison](data/selenium-reversible-net.json)

The exact precursor and product stereochemistry match the model. Neutral hydrogen selenide and acetic acid connect through explicit proton-balanced equations. The KEGG title mentions sulfide while its compound equation specifies selenide; that conflict is retained in reaction details. The selected 11-step certificate consumes three CO₂ per selenocysteine, uses permissive upstream directions and redox inputs, and permits regenerated pre-existing pools. Selenomethionine and selenohomocysteine remain unresolved.
