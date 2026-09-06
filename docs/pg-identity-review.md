# PG identity conflicts and reaction hypotheses

The main inventory remains at **2,638 / 6,220 conditional CO₂-conversion certificates**. These are balanced net-conversion hypotheses, not demonstrated Cannabis pathways.

A separate comparison examines 23 phosphatidylglycerol (PG) records whose CannabisDB names specify alkene geometry that their encoded SMILES and InChI omit. The original structures are retained. The alternatives assign only the geometry explicitly stated in the names; they are not corrections to the source records.

| Identity assertion | Reaction model | Conditional certificates |
| --- | --- | --- |
| [Original encoded structures](net.html?scenario=pg-named&comparison=original_result) | Extended chemistry | 0 / 23 |
| [Name-derived alternatives](net.html?scenario=pg-named&comparison=alternative_baseline_result) | Baseline chemistry | 0 / 23 |
| [Name-derived alternatives](net.html?scenario=pg-named) | Extended chemistry | 23 / 23 |

The extension adds 37 distinct balanced equations; 85 other proposals join existing equations. All reactants and coproducts are explicit. No identity-conversion edge is introduced. CO₂ remains the only external carbon input, with unchanged direction exclusions. Regenerated pre-existing pools are allowed; startup, energetics, transport, compartments and Cannabis enzyme activity remain unestablished.

The original 23 structures remain solver-reported infeasible in this model, not biologically impossible. Their five unresolved acyl inputs omit alkene geometry; defined alternatives are separate identities. Resolving a name/structure conflict requires source curation or analytical evidence, not an automatic stereochemistry assignment.

All 23 XML descriptions say the compounds are expected in Cannabis. That is not a reported detection or assay.

## Evidence

- [All-input closure audit](data/glycerophospholipid-input-audit.json): 21 frontier inputs, their net-supply tests and affected targets.
- [Acyl-input identity review](data/acyl-input-identity-review.json): five unresolved identities and separately encoded comparison candidates.
- [Raw XML source audit](data/pg-source-identity-audit.json): exact source fields, fragment hashes, occurrence assertions and conflicts.
- [Name-derived alternatives](data/pg-named-alternatives.json): explicit chain positions and E/Z assignments without changing the original records.
- [Reaction proposals](data/pg-named-reactions.json): exact, balanced source-template instantiations and provenance.
- [Paired net calculations](data/pg-named-net.json): every result, exact certificate and model boundary.

## Discriminating next tests

1. Resolve whether each source name or encoded structure is intended, retaining both assertions until adjudicated.
2. Establish actual Cannabis occurrence and resolve acyl positions and double-bond geometry using suitable standards and structural analysis; a generic lipid-class signal is insufficient.
3. Test the proposed substrate-specific PG/PGP synthesis and hydrolysis steps. Generic reaction templates alone do not establish activity on these exact substrates in Cannabis.

Neither positive alternative certificates nor negative model results are evidence of in vivo flux.
