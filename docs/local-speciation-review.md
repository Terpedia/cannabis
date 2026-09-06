# Local speciation sensitivity

[Explore the on-demand map, starting with glycine](net.html?scenario=speciation&target=CDB006139).

The model covers **2,723 of 6,220 historical CannabisDB records (43.8%)** with conditional, exactly replayed CO₂ net-conversion certificates. It does not establish that these pathways operate in Cannabis.

This update adds 92 balanced speciation hypotheses covering 93 previously unmatched records. Fifty equations represent net-zero intramolecular proton relocation, including neutral/zwitterionic forms. Forty-two exchange explicit H+. The earlier source-mapped protonation pass deliberately excluded the net-zero case. These new proposals are computed structural hypotheses, not curated reactions or exact ChEBI endpoint mappings.

Each pair retains heavy atoms, isotopes, bonds and encoded stereochemistry. Exact compound identities remain separate. No protein assignment is required or inferred. All 92 equations are treated as reversible sensitivity assumptions; tissue pH, pKa, compartments and physiological relevance remain unresolved.

## Reaction completeness

| Metric | Previous model | Local speciation model |
| --- | ---: | ---: |
| Historical records | 6,220 | 6,220 |
| Records with conditional CO₂ certificates | 2,638 | 2,723 |
| Distinct structural certificates | 2,636 | 2,720 |
| Records without a producing equation | 3,097 | 3,004 |
| Solver-reported infeasible records | 450 | 458 |
| Exchange-species records | 35 | 35 |
| Records participating in balanced equations | 3,115 | 3,208 |
| Balanced equations | 18,082 | 18,174 |

All earlier certificates were replayed and retained. The 85 newly covered records are represented by 84 exact structural certificates. CO₂ remains the sole external carbon input and inherited reaction-direction exclusions remain unchanged. Regenerated pre-existing pools are permitted; zero-pool startup, energetics, compartments, transport, enzyme activity and in vivo flux are not established.

All 93 targeted records gain producing equations, but eight still lack net conversion in this model: two LPA records, selenocysteine, selenomethionine, selenohomocysteine, adenylyl-molybdopterin, symmetric dimethylarginine and 1-methylhistidine. This explains the increase from 450 to 458 infeasible records; no previous certificate was lost. A producing equation alone is not a complete route.

## Evidence and next tests

- [Gap queue before this extension](data/current-reaction-gaps.json): all 3,097 prior no-producer records, exact identities and diagnostic alternatives. Charge/stereo fingerprints are search leads only.
- [Speciation hypotheses](data/local-speciation-hypotheses.json): both exact endpoint checks, explicit equations and partner-producing reaction provenance.
- [Full-inventory calculation](data/local-speciation-net.json): all records, added equations, direction exclusions and new certificates.

Review relevant compartment pH and pKa before assigning biological meaning to proton-state exchange. For the eight unresolved targets, trace every required precursor under the same carbon boundary before proposing additional reactions. Resolve source identity conflicts separately; neither a name match nor this sensitivity analysis authorizes overwriting a source structure.

## Map downloads

The target index is approximately 4.2 MB. Selecting a covered target loads shared chemistry once (14.9 MB) and its exact certificate (less than 44 KB). Gap-only selections do not need the chemistry download. Arrows retain complete input and coproduct groups; they are reaction projections, not atom flow or an execution sequence. Shared data and certificate hashes are checked before display.
