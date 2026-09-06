# Protonation review for remaining reaction gaps

The expanded structural audit covers the 17,120-equation cardiolipin model and all 6,220 target records. It identifies 557 distinct proton-transfer proposals affecting 558 records, including 381 records without a net-producing equation. These are review candidates, not additional CO₂ certificates.

## Source convention

Rhea represents reaction participants using ChEBI entities for the major microspecies at pH 7.3. Other protonation states retain separate ChEBI identities. Rhea searches can follow major-microspecies relationships even when an exact compound search is requested. Therefore a search hit alone is insufficient evidence that the queried exact structure participates in the returned reaction.

Sources checked 2026-09-06:

- [Rhea reaction participants](https://www.rhea-db.org/help/reaction-participant)
- [Rhea search and protonation relationships](https://www.rhea-db.org/help/searching-rhea)
- [Rhea downloads](https://www.rhea-db.org/help/download)

Rhea supplies [chebi_pH7_3_mapping.tsv](https://ftp.expasy.org/databases/rhea/tsv/chebi_pH7_3_mapping.tsv), containing the original ChEBI identifier, major-microspecies identifier, and whether the mapping originated from computation or curation. The [participant SMILES table](https://ftp.expasy.org/databases/rhea/tsv/rhea-chebi-smiles.tsv) validates participant structures in the source join described below.

## Exact source-join results

The Rhea tables have now been downloaded and checksummed. The participant-only table initially lacked target endpoint structures, so 393 explicit ChEBI identifier leads were retrieved from the public API. All requests succeeded; 355 records provided eligible exact-ID structures, 37 were invalid or generic for this join, and one redirected identifier remains under review.

Exact charged isomeric-structure comparison supports 302 of the 557 proposals through 305 explicit Rhea mapping rows, all marked `computation`, not `curation`. They concern 303 CannabisDB records: 289 without producing equations, nine already holding conditional net certificates, and five exchange species. Forty proposals still lack an exact target endpoint and 215 lack both endpoints in the source structures used here. Neither these identity joins nor the mapping convention establishes tissue pH, compound-specific pKa, physiological direction, or a CO₂ route.

The immutable reports are `phase1-protonation-source-join.json` and `phase1-protonation-verified-join.json`. The original pre-retrieval join remains preserved separately. No pathway metric has been increased by either join.

## Completed chemical sensitivity

The source-mapped sensitivity adds 295 distinct balanced proton-transfer equations; seven proposals join equations already present. All inherited direction exclusions and external exchanges are unchanged. The 17,415-equation model gives 259 new exact-structure certificates covering 260 additional target records, for 2,414/6,220 conditional records (38.8%). All new certificates passed independent exact net replay. The remaining inventory comprises 3,409 records without producing equations, 362 solver-reported infeasible records, and 35 exchange species. Balanced participation covers 2,803 records and is not equivalent to full net-conversion coverage.

The increase from 333 to 362 infeasible records reflects 29 previously missing-producer records acquiring equations but not feasible net conversions; no previously certified target lost its certificate. Protonation identities remain distinct, and all 305 source mapping rows retain their computed origin. This is a separate sensitivity scenario, not evidence that these equilibria or upstream pathways operate in Cannabis.

`phase1-source-mapped-protonation-net.json` is stored as 25,620 lossless rows in `terpedia-489015.terpedia_core.cannabis_phase1_source_mapped_protonation_net_20260906_v1`. Full readback matched the exported rows and SHA-256 `7a02521e69cf157e37e7c504e8c66b095098a1289b7f261f3ad7a749d843e5a4`. The static graph preserves each complete equation and its source-mapping provenance.

## Audit requirements and unresolved biological checks

1. Snapshot and checksum the mapping and structural records. Preserve the mapping's computation/curation origin.
2. Resolve both endpoints by exact structures and ChEBI identifiers; do not assign a ChEBI identifier through name-only similarity. Preserve encoded stereochemistry and isotopes.
3. Separate explicit source-mapped pairs, structure-only pairs, mismatches, and unresolved identities. A generic participant mapping does not establish exact-substrate specificity.
4. Retain balanced proton-transfer equations between distinct compound nodes. Do not merge salts, tautomers, stereoisomers, or protonation states.
5. If evaluating a chemistry-only sensitivity model, explicitly label assumed proton exchange and reversible acid–base steps. Keep parent direction exclusions and the CO₂-only external carbon boundary unchanged, and replay every new net certificate.
6. Keep tissue/compartment pH, compound-specific pKa, attainable concentrations, transport, and physiological flux unresolved unless separately evidenced. Rhea's pH convention is not evidence that every Cannabis compartment has that pH.

The structural audit deliberately excludes net-zero intramolecular proton relocation. Missing bridges are not evidence that equilibration is impossible. This review neither assigns enzymes nor establishes zero-pool pathway startup.
