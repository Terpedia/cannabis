# Defined-medium workstream

User requirement: identify non-carbon inputs and establish a minimum defined medium for the cannabis metabolic model.

## Verified starting point

`phase1-medium-inventory.json` audits every saved certificate in the selenium-forward scenario: 2,721 structures / 2,724 historical records. The full boundary contains **102 external species: CO₂ plus 101 carbon-free species**. Only **23 species (CO₂ plus 22 carbon-free species)** have positive net consumption in these selected certificates. The **35** inventory records classified as exchange species are a different metric and must not be used as the medium size.

The report retains exact charged structures, element counts, normalized input amounts, consuming targets, external outputs and provenance. It distinguishes unused and export-only species from consumed ones. No input has been proven essential or minimal by this audit.

Important boundary problems include imported iron–sulfur clusters, hydrogen peroxide and multiple phosphorus/redox forms. The selected certificates also consume N₂ in some cases; this is a model dependency to investigate, not evidence that Cannabis assimilates atmospheric nitrogen. Zero net consumption of potassium, magnesium or another ion cannot establish that the plant does not need it: the current certificates do not model biomass dilution or all catalytic requirements.

## Required model separation

1. Keep the current permissive chemistry scenario reproducible as an upper bound.
2. Curate allowed nutrient forms, gas exchange and proton/water bookkeeping separately from internal cofactors, redox carriers and reactive intermediates. Record evidence for uptake/assimilation; exact formula equality is not an identity join.
3. Retain environmental/exposure and inorganic inventory entries, but distinguish uptake from de novo synthesis. A medium that permits every detected foreign element is not automatically a normal plant growth medium.
4. Model uptake and secretion independently. Disabling an input must not inadvertently forbid disposal of that compound.
5. Optimize the shared input set for the requested set of synthesis targets, not just the intersection/union of one chosen route per metabolite. Report species-cardinality versus ingredient/salt-cardinality explicitly. Alternative optimal media may exist.
6. Distinguish a feasible set, an inclusion-minimal set (no single deletion retains feasibility), and a proven minimum-cardinality set. Greedy deletions do not prove minimum cardinality.
7. Replay exact full-input certificates under the selected boundary. Report lost coverage and missing assimilation/cofactor-regeneration reactions, rather than importing their products silently.
8. Treat light/energy, initial pools, compartments, transport, catalytic minerals, biomass demand, counterions, pH and concentrations as explicit outstanding requirements before calling it an experimentally usable plant growth medium. This workstream is not a dosing or fertilizer recipe.

The all-metabolite pathway objective remains active; minimum-medium analysis must not redefine success as supporting only today's covered subset. The covered subset can be an explicitly labeled intermediate benchmark.
