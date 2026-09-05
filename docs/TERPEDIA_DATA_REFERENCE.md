# Terpedia data reference

This file is the source contract for the cannabis carbon-provenance project.
It records what has actually been imported, how entities are identified, and
which fields are required before a reaction can contribute carbon provenance.

## Source identity warning

An accession-by-accession audit found **71 encoded-structure disagreements**
between the SDF-derived graph and CannabisDB XML; 64 occur among current
no-producer targets. Historical graphs and carbon counts describe the retained
SDF-derived structures, not independently validated name-to-structure identity.
For example, CDB006156 is named Glycerol but its SDF-derived structure is
C14H22O, while XML supplies `OCC(O)CO` and C3H8O3. XML also retains a conflicting
IUPAC name and external identifier, so wholesale replacement from XML would
not resolve every source assertion. Both assertions are preserved in
`phase1-no-producer-audit.json`; no structure or historical carbon count has
been silently replaced. Resolve these identities before interpreting their
named pathways. Passing the artifact-accounting gate does not resolve this
source-identity problem.

## Full balanced-catalog net diagnostic

The [chemistry-only net graph](net.html?scenario=catalog) keeps all 6,220
CannabisDB records searchable and admits the full 13,995 balanced-equation
catalog for a separate diagnostic. It does not change the 101-record
candidate-enzyme baseline or the 112-record inferred-completion sensitivity.

| Diagnostic result | Target records |
| --- | ---: |
| Exact net-conversion hypothesis | 304 |
| No net-producing equation in either hypothetical direction | 5,777 |
| Solver-reported infeasible under the exchange boundary | 104 |
| Explicit exchange species, not a synthesis target | 35 |

The 304 records correspond to 303 exact structures. Existing baseline
certificates are retained unchanged. The 203 additional target records require
reactions lacking baseline candidate-enzyme evidence. Across the selected
certificates, 465 distinct equations have that gap. Red directed edges identify
them; the target filter selects affected records, and highlighting dims other
edges without hiding compounds, cofactors, or coproducts. Selecting a reaction
shows the complete equation, relative extent, source links, and either its
candidate evidence or an explicit missing-evidence warning.

These are exact net stoichiometric certificates, not demonstrated Cannabis
pathways or zero-pool startup sequences. CO₂ remains the only carbon exchange;
the same 101 carbon-free exchange species and hypothetical bidirectionality
are retained. Internal pools may pre-exist but cannot be depleted overall.
Energy, thermodynamics, physiological direction, compartment compatibility,
enzyme activity and pool origin remain unverified. Inferred completion and
proton-transfer hypotheses and subsequent archive evidence supplements are
not silently added to this comparison. Atom tracing remains deferred.

The ranked gap list counts membership in selected certificates, not reaction
necessity, a minimum gene set, or a guaranteed gain from finding one enzyme.
Each gap calls for exact reaction-reference review, classification as enzymatic
versus spontaneous/catalog chemistry, correctly linked reference-protein
searches, and exact-substrate/product assays. The both-directions reactant
denominator is 11,162 distinct compounds, separate from the target denominator.

Download the [full report and ranked gaps](data/catalog-net-view/bundle.json).
Its SHA-256 is `4fdc58d257093a70c9b495b865d5b97f1477ed66b4a49887abcf89f903860e14`.
The static bundle is byte-identical to `data/reports/phase1-catalog-net-gaps.json`
and pins seven source reports. Independent tests replay all positive certificates
with rational arithmetic, audit full-equation isotope/element/charge balance,
verify target identity, candidate evidence and gap membership, and check
no-producer statuses against the full catalog. Infeasibility remains a retained
numerical solver result, not an exact impossibility proof.

Reproduce with `PYTHONPATH=src python -m cannabis_carbon.phase1_catalog_net_gaps`
then `PYTHONPATH=src python -m cannabis_carbon.phase1_catalog_net_view`.
The full-catalog solve may take several minutes. The static report is stored
only once in Pages to avoid duplicate large downloads.

GCP table `terpedia-489015.terpedia_core.cannabis_phase1_catalog_net_gaps_20260904_v1`
contains 9,159 verified records: 6,220 targets, 303 certificates, 880 selected
reactions, 855 compounds, 465 ranked gaps, 435 evidence records and one metadata
row. The destination was absent before loading; every complete stored record
was read back and matched exactly to its local export. The service-account
identity, load job and verification are recorded in
`data/reports/phase1-catalog-net-gaps-gcp.json`.

## Protein discovery for previously unscreened catalog gaps

The 465 selected-certificate enzyme gaps were compared by exact reaction ID
against both earlier full-proteome screens. Existing results for 266 equations
are retained; 199 equations were not previously screened. No failed or weak
search was silently promoted or resubmitted. The new queue preserves full
equations, original source IDs, and selected-target membership.

Eight successful reviewed, non-fragment UniProt lookups using explicit
published Rhea families returned 2,184 reference proteins for 166 equations;
33 equations returned no reviewed reference. All 2,184 sequences were retrieved
and searched against all 30,304 Cannabis reference-proteome proteins with
DIAMOND sensitive mode. Thresholds remain E-value ≤ 1e-5, identity ≥ 30%, and
both query and reference coverage ≥ 50%; no target-count truncation was used.

Of 9,452 raw alignments, 2,924 passed. They yield 296 Cannabis proteins and
595 protein–reaction hypotheses across 97 equations. The other new equations
retain 44 weak-hit-only, 25 no-hit and 33 no-reference-sequence results.
The reference records, full sequences, alignments and retrieval snapshots
retain checksums and reaction-level joins. Candidate status is not an assay,
exact substrate specificity, physiological direction, or tissue compatibility.

As a separate evidence comparison, these 97 equations occur in selected net
certificates for 181 target records. Of the 203 records with enzyme gaps,
only **uric acid (CDB004839)** has all its selected missing steps covered by
the new candidates. Across the catalog diagnostic, 368 distinct selected gap
equations still lack candidates. This comparison does not establish a new
physiological pathway, zero-pool startup, or necessity of a selected step.
The baseline and original catalog chemistry snapshots remain unchanged. The
first catalog evidence supplement added candidates to 97 reactions, leaving
368 gaps. The current combined supplement below includes the backfill results;
the affected-target filter uses the remaining gaps. All full equations,
extents, participants, external exchanges and startup results are unchanged.
Atom tracing remains deferred.

The [earlier evidence supplement](data/catalog-net-view/evidence-v1.json)
has SHA-256 `ff242379d6c77127f112a8e170b4210284ef78b7e94f174d1faf30860dff0840`.
It contains 97 candidate-evidence records, 181 exact-structure certificate
updates and 181 target updates, plus metadata. It retains the original missing
reaction lists alongside the remaining lists. The net graph now reports 102
selected target certificates with candidates for all steps and 202 with gaps,
out of the unchanged 304 net-conversion hypotheses. This is evidence coverage,
not established pathway completeness. A failed or mismatched supplement load
stops rendering and offers a retry rather than presenting mixed snapshots.

Reproduce the supplement with
`PYTHONPATH=src python -m cannabis_carbon.phase1_catalog_evidence`, followed by
`PYTHONPATH=src python -m cannabis_carbon.phase1_catalog_net_view`. Pages retains
the original catalog bundle once and fetches the 2.6 MB supplement separately;
it does not duplicate the full net graph. The manifest pins both reports.
GCP table `terpedia-489015.terpedia_core.cannabis_phase1_catalog_evidence_20260904_v1`
contains all 460 supplement records, verified by complete-record readback.
The receipt is `data/reports/phase1-catalog-evidence-gcp.json`.

Reproduce reference discovery with
`PYTHONPATH=src python -m cannabis_carbon.phase1_catalog_references` and screening
with `PYTHONPATH=src python -m cannabis_carbon.phase1_catalog_protein_search`.
Reference lookup caches retain unsuccessful retrievals explicitly. Scientific
tests replay exact-family joins, the prior/new queue partition, sequence hashes,
alignment thresholds, reaction joins and raw-hit outcomes.

| Report | SHA-256 | Verified GCP rows |
| --- | --- | ---: |
| `phase1-catalog-references.json` | `cd230249677b4c5775d910777f285be2218c9d0c91ba6177908d7c33120f97e9` | 2,857 |
| `phase1-catalog-protein-search.json` | `2f31678e51e644a5787cfb4fdf8cd14131f837ae06d139575cff20b5706025e9` | 5,641 |

Tables in `terpedia-489015.terpedia_core` are
`cannabis_phase1_catalog_references_20260904_v1` and
`cannabis_phase1_catalog_protein_search_20260904_v1`. Every complete stored
record was read back and matched exactly to the local export. Load jobs and
service-account verification are recorded in
`data/reports/phase1-catalog-protein-gcp.json`.

## Backfill of previously skipped reference-family lookups

The 368 reaction gaps remaining after the catalog supplement were audited
against all three protein screens and their actual lookup histories. Of these,
211 had reference sequences already screened. The other 157 split into **51
families never queried**, 98 with completed reviewed-family lookups, and eight
without a published family mapping in the pinned Rhea snapshot.
This distinguishes `no-reference-sequence` from a completed negative lookup;
the earlier priority pass had not queried every reaction family.

Three new reviewed, non-fragment lookup batches covered exactly those 51
unqueried master families. Existing successful lookups were reused with
checksum verification. The combined reference inventory yields 431 proteins
for 43 equations; 106 still have no reviewed reference returned, and eight
still have no published family mapping. No unreviewed annotations or numeric
Rhea-ID arithmetic were introduced.

All 431 reference sequences were retrieved and screened against all 30,304
Cannabis reference-proteome proteins using the unchanged DIAMOND thresholds.
Of 1,553 raw alignments, 375 passed, yielding 111 candidate proteins and 132
protein–reaction hypotheses for **19 reactions**. Outcomes for the other
backfill records are 12 weak-hit-only, 12 no-hit and 114 no-reference-sequence.
These are candidate leads, not demonstrated catalytic activity, specificity,
physiological direction or complete Cannabis pathways.

An evidence-only comparison reduces the distinct remaining reaction gaps
from 368 to 349, but **does not close any additional selected net certificate**.
The catalog graph now incorporates this backfill through the separately
versioned combined supplement below; its earlier 368-gap evidence snapshot
remains downloadable. Original chemistry, atom
accounting and all prior search results are preserved; atom tracing remains
deferred.

Reproduce using `PYTHONPATH=src python -m cannabis_carbon.phase1_reference_backfill`
and `PYTHONPATH=src python -m cannabis_carbon.phase1_backfill_protein_search`.
The reports preserve all 368 gap classifications, exact source/family joins,
lookup requests and cached responses, full sequences, and checked alignment
results. Required snapshots are committed for clean-checkout replay.

| Report | SHA-256 | Export rows |
| --- | --- | ---: |
| `phase1-reference-backfill.json` | `c135d4f5563e7b9d0d36c80d4227e6aa0898ee7e7156207f9634c622d3c99508` | 1,004 |
| `phase1-backfill-protein-search.json` | `c47ad4f250762b2f835827260adf403ea032c540bc1ca75f24e8f5dcdb4ff0ef` | 1,083 |

Both exports are verified in `terpedia-489015.terpedia_core`, in tables
`cannabis_phase1_reference_backfill_20260904_v1` and
`cannabis_phase1_backfill_protein_search_20260904_v1`. Every complete record
was read back and matched to its local export. The receipt is
`data/reports/phase1-reference-backfill-gcp.json`.

### Combined catalog evidence snapshot

The [current graph evidence supplement](data/catalog-net-view/evidence.json)
preserves all 97 earlier reaction-evidence records verbatim and adds the 19
backfill candidates. Its 116 reaction records contain 385 distinct Cannabis
candidate proteins (deduplicated across screens). There are 349 remaining
reaction gaps. Candidate evidence improves for 191 target records relative to
the original catalog diagnostic, but the number of selected certificates with
candidates for every step stays at 102; 202 still have gaps.

This is a union applied once to the original catalog bundle, not a second
application on top of the first supplement. All equations, coefficients,
relative extents, metabolites, exchange assumptions, and startup labels remain
unchanged. Each backfill record retains its own search-report link and evidence
class. Blue arrows denote screened candidates, not demonstrated activity;
red arrows identify the remaining enzyme-candidate gaps.

SHA-256: `9cf7e4c5e9b8b47058aeb90eded2f5bd8faa4741901178630c43a8b08860435e`.
The report `data/reports/phase1-combined-catalog-evidence.json` contains 499
export records: 116 enzyme-evidence records, 191 certificate updates, 191
target updates and one metadata record. The previous report remains unchanged
at `data/reports/phase1-catalog-evidence.json` and downloadable as
`data/catalog-net-view/evidence-v1.json`.

Reproduce with `PYTHONPATH=src python -m cannabis_carbon.phase1_combined_catalog_evidence`
then `PYTHONPATH=src python -m cannabis_carbon.phase1_catalog_net_view`.
Tests check exact prior-record preservation, duplicate rejection, complete gap
recalculation, source hashes, static data identity and all graph participants.
GCP table `terpedia-489015.terpedia_core.cannabis_phase1_combined_catalog_evidence_20260904_v1`
contains all 499 records, verified by full-record readback. The receipt is
`data/reports/phase1-combined-catalog-evidence-gcp.json`.

## Imported Terpedia snapshot

Source repository: [Terpedia/terpedia-knowledge](https://github.com/Terpedia/terpedia-knowledge)

Organism: *Cannabis sativa*

Reference proteome: `https://rest.uniprot.org/proteomes/UP000583929`

Snapshot files in this repository:

- `data/terpedia/cannabis-sativa-metabolic-network.json.gz`
- `data/terpedia/cannabis-sativa-metabolic-catalog.json`
- `data/terpedia/cannabis-sativa-metabolic-reachability.json`
- `data/terpedia/cannabis-sativa-metabolic-hypotheses.json`

Snapshot inventory:

| Entity or relationship | Count |
| --- | ---: |
| proteins | 30,304 |
| metabolites | 1,245 |
| biochemical reactions | 1,086 |
| enzyme classes | 765 |
| `has_reactant` statements | 2,353 |
| `has_product` statements | 2,783 |
| `catalyzes` statements | 4,337 |
| `annotated_with_ec` statements | 3,630 |

The snapshot is a versioned input, not a claim that the network is complete or
that every annotation represents activity in cannabis tissue.

The working network also includes the source-linked curation addition
`data/terpedia/reaction-additions.json`: Rhea 13429 / directional Rhea 13430
(IspD, EC 2.7.7.60), connecting MEP-4P to CDP-ME. This adds reaction knowledge
for network completeness and balance auditing; it is not evidence of a
cannabis IspD enzyme or in-vivo flux.
The same curation file adds three exact-structure cannabinoid oxidocyclization
reactions: CBGA to CBDA, THCA-A, and CBCA, each with O₂ and H₂O₂ explicitly
represented. For CBDA synthase, the substrate/product conversion is directly
characterized but the O₂/H₂O₂ mechanism has conflicting historical evidence;
that uncertainty is retained on the reaction record. These are curation edges,
not proof that the corresponding cannabis proteins are present or active in
every accession.
The overlay `data/terpedia/varin-reaction-additions.json` adds CannabisDB structures
for CBGVA, CBGV, CBCVA, and CBCV, plus balanced candidate edges for the varin
prenylation/oxidocyclization branch and the corresponding acid decarboxylations.
It also adds the balanced coupled TKS/OAC hypothesis from butyryl-CoA and three
malonyl-CoA molecules to divarinolic acid, 4 CoA, and 3 CO₂, with separate TKS
and OAC genome candidates. The reaction is a coupled bookkeeping edge, not a
claim of a single multifunctional enzyme.
CBGVA-to-CBDVA and CBGVA-to-THCVA use EC-number joins to surface genome candidates;
the edges remain hypotheses because the cited yeast reconstruction does not prove
the exact Cannabis in-planta enzyme route.
It also adds the exact ChEBI:66955 olivetolic acid node and the
geranyl-diphosphate transfer reaction that produces CBGA; the olivetolic acid
structure is sourced from ChEBI/PubChem and the reaction from the hemp
geranylpyrophosphate:olivetolate geranyltransferase report.
It also adds the upstream type-III polyketide step from hexanoyl-CoA and three
malonyl-CoA to 3,5,7-trioxododecanoyl-CoA plus CO₂, and the exact Rhea 34123
OAC cyclization to olivetolate. Protonation is explicit so both additions are
auditable for atoms and charge.

## Entity identifiers

- Metabolites: ChEBI IDs such as `chebi:15377`.
- Reactions: Rhea master IDs such as `rhea:25500`.
- Proteins: UniProt IDs such as `uniprot:a0a1v0qsg0`.
- Enzyme classes: EC identifiers represented by the network's entity record.

Never merge entities by label alone. Preserve exact stereochemistry, charge,
protonation, salts, adducts, and organism/source identifiers.

## Reaction fields

Each reaction entity may contain:

- `attributes.equation`: human-readable biochemical equation.
- `attributes.reactionSmiles`: structure-bearing reaction input for RDKit.
- `attributes.reactionSmilesRheaId`: directional Rhea record when present.
- `attributes.ecNumbers`: exact EC numbers associated with the reaction.
- `elementBalance` and `chargeBalance`: source audit fields.
- `url`: source record URL.

Participant statements use:

- `subjectId`: reaction ID.
- `predicate`: `has_reactant` or `has_product`.
- `objectEntityId`: exact metabolite ID.
- `qualifiers.stoichiometricCoefficient`.
- `qualifiers.compartment`.

Direction predicates are evidence-bearing, not assumptions:

- `physiological_direction_left_to_right`
- `physiological_direction_right_to_left`

## Enzyme evidence

An enzyme association is represented by `catalyzes`, `maps_to_reaction`, or
Terpedia's source-native `has_catalytic_activity` predicate. The exporters
normalize all three forms to reaction enzyme associations.
The associated statement's qualifiers must be retained, especially:

- `directExperimentalEvidence`
- `inference`
- `warning`
- `assertionStatus`
- `evidenceBoundary`
- `sources`

An EC-to-Rhea join or genome homology hit is a candidate producer, not a
confirmed cannabis enzyme. Candidate-protein discovery must attach the genome
assembly, protein ID, search method, score, alignment coverage, motif/domain
evidence, and source URLs.

## Carbon provenance contract

For every product carbon atom, the mapping output must contain either:

1. a precursor carbon atom and the reaction/provenance edge that carries it;
2. a candidate precursor mapping with status `candidate` or `inferred`; or
3. an `unresolved` record with the reaction ID, product atom, and blocking cause.

Structural RDKit mapping is not isotope tracing and does not establish flux.
New carbon introduced by fixation, transfer, methylation, carboxylation, or an
unmodeled co-substrate must not be assigned to an arbitrary reactant carbon.

## Current completeness baseline

Atom tracing is deferred for Phase 1. Acceptance requires full element and charge
balance; candidate enzyme evidence and CO₂ connectivity are separate metrics.

### Original-source archive recovery

The [archive resolution report](data/phase1-archived-references.json) recovers all
47 requested original-source sequences: 46 UniParc identifiers in Terpedia's
MARTS records plus the explicit archive link for inactive UniProt A0A4P2VJ76.
The latter is accepted only after the archived record independently retains
the old accession as a cross-reference. Its inactive status is not rewritten.
All 104 source associations across 66 completion equations are retained,
including four CannabisDB target records. Raw responses, retrieval dates,
SHA-256 hashes, sequence lengths, MD5 checks and complete archive cross-references
are preserved. No active accession or functional annotation is substituted.
[UniParc is a sequence archive](https://www.uniprot.org/uniparc), not evidence of
exact enzymatic activity.

The [archived-reference screen](data/phase1-archived-protein-search.json) searches
all 30,304 Cannabis proteins and retains 107 raw alignments, 65 passing alignments,
44 distinct Cannabis candidates and 65 protein–reaction hypotheses across three
equations. The same identity (≥30%), two-sided coverage (≥50% each) and E-value
(≤1e-5) thresholds apply. Of the 66 equations searched, 57 have no hits and six
have only weak hits; these outcomes are not evidence of biological absence.

The [evidence supplement](data/phase1-archived-evidence.json) adds first candidate
leads for two completion equations and additional sequence support for one that
already had leads. One newly linked equation involves CannabisDB CDB004889,
(S)-2,3-epoxysqualene. Participation in that equation is not a demonstrated
producing pathway: exact identities, full inputs, direction and original source
product warnings remain explicit. No previously unsupported priority target
receives a passing archived-reference candidate in this screen.

The [completion graph](completions.html) shows the archive screen separately
inside protein evidence and includes its passing leads in hypothesis filters.
Across all 765 completion equations, candidate-lead coverage becomes 387
(previously 385); 378 remain without a lead. The 62 target records with at least
one lead and five targets with no lead remain unchanged. Original chemistry and
protein reports are immutable. **CO₂ baseline and completion-sensitivity
certificates are not changed or re-scored by this evidence-only supplement.**
Atom tracing remains deferred.

Reproduce in order with `PYTHONPATH=src python -m cannabis_carbon.` plus
`phase1_archived_references`, `phase1_archived_protein_search`,
`phase1_archived_evidence`, and `phase1_completion_view` (join the prefix and
module name into one argument). The reference resolver reuses checksummed cached
responses; it does not silently replace a failed or archived record.

| Report | SHA-256 | Prepared Terpedia rows |
| --- | --- | ---: |
| `phase1-archived-references.json` | `eb2d4a7e538dbd063fec07252bcff7be972e21a3e07ff8f4b50f7ba149049ae5` | 161 |
| `phase1-archived-protein-search.json` | `c77f50fd2b85009c6ca103a3dc8960aece7971cfcd4cea87cbac6dc35ece5399` | 270 |
| `phase1-archived-evidence.json` | `6da25f45823409f07c7ca2822e33a7cef58be08d3c05a1475e681f22ac4cd6e0` | 67 |

These three exports were loaded and fully read-back verified on 2026-09-04
using `cannabis-metabolome@terpedia-489015.iam.gserviceaccount.com`.
Tables in `terpedia-489015.terpedia_core` are
`cannabis_phase1_archived_references_20260904_v1`,
`cannabis_phase1_archived_protein_search_20260904_v1`, and
`cannabis_phase1_archived_evidence_20260904_v1`. Each destination was absent
before loading. All 498 stored records match the local exports exactly,
including complete JSON payloads and report hashes. Job IDs and verification
details are retained in `data/reports/phase1-archived-gcp.json`.

### Completion-network sensitivity: 11 additional net-balance hypotheses

The [separate sensitivity graph](net.html?scenario=completions) augments the
1,472 baseline candidate equations with 321 original-MARTS-homology completion
hypotheses. Of these additions, 146 were already balanced equations in the full
Terpedia catalog but lacked baseline candidate support, and 175 are additional
inferred equations. The 64 completions already carrying baseline candidate
evidence are not double-counted. The 380 completions without a candidate lead
remain excluded. All 6,220 target records are retained.

With the **same exact CO₂ and carbon-free external exchanges**, 11 additional
target structures obtain rationally reconstructed net-balance certificates:
112 target records / 111 unique structures versus the unchanged baseline of
101 records / 100 structures. Each new certificate has 22–23 directed reaction
steps and depends on one newly admitted completion. Every coefficient and
coproduct is retained, every internal species has nonnegative net production,
and CO₂ supplies all net carbon consumed. Both directions remain hypothetical;
positive results are conditional on unverified chemistry, enzyme specificity,
compartments, energy and thermodynamics, not confirmed Cannabis pathways.

| Newly conditionally feasible target | CannabisDB ID |
| --- | --- |
| Farnesol | CDB006149 |
| alpha-Cubebene | CDB000109 |
| alpha-Bisabolol | CDB000119 |
| alpha-Phellandrene | CDB000199 |
| beta-Phellandrene | CDB000200 |
| alpha-Thujene | CDB000206 |
| trans-Sabinene hydrate | CDB000302 |
| Ledol | CDB000316 |
| Gamma-Muurolene | CDB000471 |
| Beta-maaliene | CDB000572 |
| Bulnesol | CDB006348 |

**Zero-pool startup does not improve** in either the CO₂-only or permissive
carbon-free scenario. The new net certificates allow regenerated, pre-existing
internal pools; they do not establish pool synthesis or origin. In the net
scenario, 187 target records remain solver-reported infeasible, 5,886 lack a
net-producing candidate equation, and 35 are explicit exchange species rather
than synthesis targets. Numerical infeasibility is not biological absence.

The baseline report and its certificates are unchanged. Orange arrows identify
completion-sensitivity edges; the graph retains complete required-input and
coproduct lists, relative extents, exact structures, homology alignments and
original-source product/stereo warnings. The graph is not a reaction execution
sequence or atom-flow map. The selected certificates are neither unique nor a
functional ranking of the proteins.

Each new certificate includes test proposals: resolve exact source-product
identity first; express candidate proteins and assay the exact original organic
substrate while identifying products and quantifying inferred inorganic
coproducts; then investigate full-route input supply, expression, compartments
and physiological direction. A homology match is not a positive assay.

Download [the full sensitivity report](data/phase1-completion-connectivity.json).
SHA-256: `79f6349b0318f0ac727fde268e34588592421ad4c470cb80e8dc4f659b554a08`.
Reproduce with `PYTHONPATH=src python -m cannabis_carbon.phase1_completion_connectivity`
and `PYTHONPATH=src python -m cannabis_carbon.phase1_completion_net_view`.
The 7,396-row Terpedia export is verified in GCP (6,220 targets, 321 admitted reactions,
380 excluded completions, 11 certificates, 24 baseline certificate reactions,
437 compounds, two startup scenarios and one metadata row). Versioned table:
`terpedia-489015.terpedia_core.cannabis_phase1_completion_connectivity_20260904_v1`.
The destination was verified absent before loading on 2026-09-04. The dedicated
Cannabis service account completed the load, and all 7,396 complete records
were read back and matched exactly to the local export. The load job and
verification are recorded in `data/reports/phase1-completion-connectivity-gcp.json`.

### Completion protein evidence

The completion graph now includes a separate, evidence-only protein overlay.
An exact full-equation join retains existing candidate evidence for 64 of 765
completion hypotheses. The remaining 701 equations were screened using only
their original MARTS protein references, never enzymes borrowed from an inorganic
stoichiometry template. Of 915 requested UniProt identifiers, 869 sequences were
retrieved and 46 unreturned identifiers remain explicit. All 30,304 proteins in
the Cannabis reference proteome were searched. Passing thresholds were identity
at least 30%, query and reference coverage each at least 50%, and E-value at most
1e-5. The screen retained 10,155 passing alignments involving 55 Cannabis proteins,
321 equations and 6,072 protein–reaction hypotheses. These are homology leads,
not verified exact-product activity or validated inorganic stoichiometry.

Across all completions, 64 have existing exact-equation candidate evidence,
321 have new original-source homology leads, and 380 have no screened candidate
lead. Of the 67 target records with completions, 62 have at least one lead and
five have none. Eight targets have an unsupported alternative, including three
mixed-evidence targets; the graph's protein filter operates on individual
hypotheses so these alternatives remain visible. Nine of the ten priority target
records have passing new candidates. No CO₂-pathway completeness metric changes.
Inactive, missing and unsupported source references and previous product/stereo
warnings remain explicit. Representative alignments are chosen by bit score for
display, not as a ranking of functional specificity.

Reports: [discovery queue](data/phase1-completion-protein-discovery.json),
[full proteome screen](data/phase1-completion-protein-search.json), and
[evidence overlay](data/phase1-completion-protein-evidence.json).
Reproduce with `PYTHONPATH=src python -m cannabis_carbon.` followed by each module
name in order: `phase1_completion_protein_discovery`,
`phase1_completion_protein_search`, `phase1_completion_protein_evidence`, and
`phase1_completion_view` (the module prefix and name form one argument).

Versioned tables in `terpedia-489015.terpedia_core`:

| Table | Rows | Report SHA-256 |
| --- | ---: | --- |
| `cannabis_phase1_completion_protein_discovery_20260904_v1` | 1,135 | `cad0025e3031ca15425d691786163cbfb419657dc6b463f38ed834a7c2db7e49` |
| `cannabis_phase1_completion_protein_search_20260904_v1` | 11,797 | `c3b1412f9901e808549db69fd9ace9e5f6d9659d41a4aa2764247d294c60197d` |
| `cannabis_phase1_completion_protein_evidence_20260904_v1` | 766 | `8b3180990d046ea52fedbe3527a6ffb4bc80d02d28e6c8f44efd149cccba815a` |

The [stoichiometric completion graph](completions.html) is a separate,
review-only Cytoscape view of [765 balanced completion hypotheses](data/phase1-marts-completions.json).
They cover 67 CannabisDB target records, including all ten newly identified
MARTS target gaps. Of the 765 distinct balanced equations, 310 already exist
in the full-network baseline and 455 are additional inferred equations. These
counts do **not** change candidate-enzyme or CO₂-pathway completeness.

The generator retains all 4,315 unbalanced source rows grouped into 834 recorded
equation variants: 765 receive a hypothesis, 17 fall outside the single-organic-
substrate/product scope, and 52 have no compatible reference template. It scans
6,114 eligible orientations of independently balanced Rhea-backed equations.
Every hypothesis preserves the original MARTS organic substrate and product,
including charge, isotope and stereochemical encoding. It copies complete
carbon-free participant lists only when the reference has the **same exact
organic substrate** and a sole organic product with the same elemental/isotopic
composition and formal charge. The reference product can be a different isomer;
that compatibility is explicitly **not an identity, mechanism or enzyme join**.
No carbon-bearing cofactor, donor or seed is inserted. Every resulting equation
is independently checked for element, isotope and charge balance.

These are inspectable hypotheses, not silent repairs: original MARTS equations
remain unbalanced in their immutable audit. Reference directions are hypothetical,
all-input supply is unestablished, and enzyme evidence remains empty in the
immutable chemistry report. Candidate leads are attached only through the
separate protein-evidence overlay described above.
Source product identity, full coproduct stoichiometry, protein specificity and
compartment conditions require review or experiments. The earlier stereochemical
and product-label warnings remain in force.

The static view retains all 6,220 target records, filters to targets with
hypotheses or all targets, and shows one complete equation at a time. Amber nodes
are inferred inorganic species; green nodes are unchanged organic compounds.
Directed dashed arrows project input/output pairs and are not additional
reactions or atom-flow claims. Full coefficients, source references and review
requirements remain available. A versioned 5.51 MB bundle contains only the
target-linked subset; the complete 16.38 MB report preserves all 765 hypotheses,
355 reference equations and 779 participating/reference compound structures.

Reproduce with `PYTHONPATH=src python -m cannabis_carbon.phase1_marts_completions`
and `PYTHONPATH=src python -m cannabis_carbon.phase1_completion_view`.
Terpedia table:
`terpedia-489015.terpedia_core.cannabis_phase1_marts_completions_20260904_v1`,
8,954 rows (834 variants, 765 completions, 6,220 targets, 355 reference reactions,
779 compounds, one metadata record). Report SHA-256:
`ecd10d6bccac52643be2eead4e6d57291d87d05ab546b4fb4d80ae698c2222fe`.
Atom tracing remains deferred.

The [MARTS gap reference review](data/phase1-marts-gap-references.json) follows
all **32 source rows for the ten newly identified stoichiometry-gap targets**.
It resolves 23 active UniProtKB records, preserves one inactive entry and its
deletion/UniParc metadata, and flags three source UniParc IDs for separate
resolution. Two MARTS rows lack a UniProt reference; their other source fields,
including any GenBank accession, are retained. The 27 requested identifiers
are not 27 characterized Cannabis enzymes.

Ten source rows have catalytic Rhea annotations. The review checks their
explicit published direction families against 36 retained Terpedia Rhea source
equations, but finds **zero balanced alternatives with exactly matching
carbon-containing participants and coefficients**. This is a scoped annotation
result, not evidence that the reactions or enzymes do not exist. All ten target
gaps remain unresolved; no enzyme evidence is transferred and no equation is
silently repaired. Full source responses, sequences, citations, evidence codes,
retrieval timestamps and checksums are preserved for further review.

Two source rows for CDB000302 (trans-sabinene hydrate), referencing
[L0HB77](https://www.uniprot.org/uniprotkb/L0HB77/entry) and
[L0HAM7](https://www.uniprot.org/uniprotkb/L0HAM7/entry), yield a **stereo-only
diagnostic lead**. The retained Rhea 19566/19567 product has unspecified
stereochemistry, whereas the CannabisDB/MARTS product encodes stereochemistry.
The diagnostic removes stereo only for comparison; it preserves isotopes,
charge, bonds and coefficients and never becomes an exact-identity or enzyme
join. Review the cited study (PMID 23246843), including relative versus absolute
product stereochemistry, before specializing a generic equation.

Another priority is CDB000078: MARTS:2138 links to
[A0A348B782 / PpSTS-06](https://www.uniprot.org/uniprotkb/A0A348B782/entry), whose
function annotation names 9-epi-caryophyllene. The cited
[primary study's Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC6116744/)
distinguishes that product from caryophyllene produced by PpSTS14. This flags a
product-label/stereochemical review; it is not an automatic correction to
CannabisDB or a claim that its encoded identity is wrong.

Reproduce with `PYTHONPATH=src python -m cannabis_carbon.phase1_marts_gap_references`.
The verified Terpedia export is
`terpedia-489015.terpedia_core.cannabis_phase1_marts_gap_references_20260904_v1`:
106 rows (10 targets, 32 source reviews, 27 reference outcomes, 36 Rhea sources,
one metadata record). Report SHA-256:
`04b948fdb153ea1d18a89afd21e0440b061361f9094debca9cf842d10e9c6328`.
Atom tracing remains deferred; all-input supply and Cannabis pathway claims
remain unestablished.

The [whole-MARTS audit](data/phase1-marts-audit.json) scans all 4,639 MARTS
records in Terpedia's normalized reaction catalog. **324 source rows balance**,
deduplicating to 73 exact full equations, including **27 absent from the frozen
full balanced-network baseline**. The report retains these as candidate additions;
it does not rewrite the existing network or its pathway certificates. Six target
records participate in balanced MARTS equations; all six already had balanced
network participation. Therefore this audit adds **zero newly covered targets**.

The other **4,315 source rows are imbalanced**, representing 834 distinct recorded
reaction strings. All source fields, protein/EC references, source orientation,
and element/charge deficits are preserved. There are 68 targets with exact
participation in these excluded equations, including **10 without baseline
balanced participation**: farnesol, beta-caryophyllene, alpha-cubebene,
alpha-bisabolol, alpha-phellandrene, beta-phellandrene, alpha-thujene,
trans-sabinene hydrate, ledol and bulnesol. These are exact encoded-structure
leads, not name-only matches or completed pathways.

For example, MARTS:1 records a diphosphate-bearing substrate turning into a
hydrocarbon but omits coproducts. The recorded right-minus-left deficit is
H −1, O −7, P −2, with charge +3. The audit **does not automatically insert
pyrophosphate, water or protons**; full stoichiometry and the source publication
must be reviewed before an equation enters a balanced pathway scenario.
The source's `directly_characterized_MARTS_record` label is retained as a source
claim, not converted into Cannabis enzyme evidence. Every new equation has
explicit full sides and coefficients, unresolved physiological direction and
all-input supply, and no automatically assigned Cannabis enzyme.

Reproduce with `PYTHONPATH=src python -m cannabis_carbon.phase1_marts_audit`.
Use `--fetch` only for the first immutable snapshot download and configure
`CANNABIS_BQ` when `bq` is not on PATH. The report pins its full-network and raw
MARTS inputs by SHA-256; the source ledger preserves every fetched row even when
the local raw cache is absent. The Terpedia export is
`terpedia-489015.terpedia_core.cannabis_phase1_marts_audit_20260904_v1`:
11,013 rows (6,220 targets, 4,639 sources, 73 reactions, 80 compounds, one
metadata record). Report SHA-256:
`1d021f0813989abcb7d7065ed8791b6013c49c2263fa6315c418d591613b94ea`.
Atom tracing remains deferred. Next: review complete equations for the ten
new target leads, then validate references and screen the Cannabis proteome.

The [protonation identity audit](data/phase1-protonation-audit.json) retains all
6,220 CannabisDB target records and proposes **342 separate, balanced
proton-transfer hypotheses for 343 target records**. Of those targets, 341
currently have no exact balanced-reaction participation. These are review leads,
not added candidate-enzyme evidence or established CO₂ pathways; all existing
coverage and net-conversion metrics remain unchanged.

For example, CDB006145 (methylamine) yields the balanced hypothesis
`CN + [H+] → C[NH3+]`. H⁺ is an explicit required input in that display
orientation; its compartment availability and the physiological direction are
unresolved. No protein, EC assignment or metabolite-observation evidence is
created by this structure comparison.

Each hypothesis retains both exact compound IDs, full reaction SMILES including
every H⁺, independent element/charge checks, source reaction participation,
normalization checks, and explicit review requirements. The audit uses RDKit's
[protonation-only Uncharger](https://www.rdkit.org/docs/cppapi/classRDKit_1_1MolStandardize_1_1Uncharger.html)
solely for candidate lookup. It independently rejects changed heavy-atom
connectivity, bond order, isotopes or stereochemistry; changed protonation sites
are restricted to N/O/P/S. It does not merge compounds, normalize tautomers,
strip salts, or infer enzyme activity. Net-zero proton relocations, retained
explicit hydrogen atoms, radicals, generic and multicomponent structures are
outside this conservative audit. There are 21 excluded target records; a
negative result means no match within these restrictions, not chemical absence.

Before a bridge enters any separately labeled pathway sensitivity scenario,
review source identity, pKa, compartment pH, direction, and whether spontaneous
equilibration or catalysis applies. CannabisDB labels are source labels: a
concrete encoded structure does not establish that a polymer or mixture label
has been resolved correctly. Atom tracing remains deferred.

Reproduce with `pip install -e '.[protonation]'` and
`PYTHONPATH=src python -m cannabis_carbon.phase1_protonation_audit`.
The report pins its full balanced-network input by SHA-256 and records RDKit
2026.3.6. Its Terpedia export is
`terpedia-489015.terpedia_core.cannabis_phase1_protonation_audit_20260904_v1`:
7,247 rows (6,220 targets, 342 bridges, 684 compounds, one metadata record).
Report SHA-256:
`7698099beebe91f63147b4f44e761b42bba5157189b7064a86a0ae6c983cedef`.

The [multi-reaction net-conversion explorer](net.html) now renders the complete
selected net-conversion certificate in static Cytoscape, or one selected reaction.
It retains all 6,220 target records, including gaps, and offers a certificate-only
filter for the 101 target records with net solutions. Default selection is a
certified limonene record; `?target=CDB…` can select another exact record.
Search and target selection preserve CannabisDB IDs instead of merging labels.

Compounds are nodes. Directed arrows project each full reaction's input/output
pairs; they are **not additional reactions, atom-flow claims, or an execution
order**. The view reports compound, directed-reaction and projection counts
separately. Each edge retains the full required inputs, outputs, coefficients,
relative extent, reaction ID, evidence IDs and candidate protein identifiers.
Selecting an edge or reaction exposes its full equation, source links and
unaltered evidence records. Zero-net internal participants are distinguished
from net inputs and products and can be highlighted without hiding any species.
The full certificate balance remains visible even while focusing on one step.

The [graph manifest](data/net-view/index.json) pins the net report and all three
source evidence catalogs by checksum. Its ~6.7 MB uncompressed bundle contains
229 referenced evidence records, exact equations, certificates and target
statuses. It is generated with
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_net_view`.
Client requests revalidate the manifest and version the bundle URL by checksum.
Automated checks cover every projected certificate, source fidelity, no omitted
participants, reaction focus, pool highlighting, gaps, retry behavior and manifest
path restrictions. No browser visual QA was requested for this release.
The underlying scientific report remains the verified GCP net-flux snapshot;
this display layer adds no new reaction, enzyme, startup or physiological claim.
Atom tracing remains deferred.

The [candidate-linked net-conversion audit](data/phase1-candidate-net-flux.json)
tests a different question from zero-pool startup: can candidate-linked balanced
equations make a target with CO₂ as the only net carbon input while fully
regenerating internal compounds? The model uses the same 1,472 candidate-linked
equations and permits CO₂ plus 101 carbon-free exchange species. Both directions
remain hypothetical. Every nonexchange compound must have nonnegative net
production; target production must be at least one unit. Every positive net
product is explicitly recorded as an export, and every net-consumed species is
listed. No ATP, NAD(P)H, CoA or other organic compound is allowed as a net input.

The results retain all 6,220 CannabisDB records:

- 101 target records (100 exact structures) have exact net-conversion certificates.
- 187 have solver-reported infeasibility under these constraints.
- 5,897 have no net-producing candidate-linked equation.
- 35 are explicit exchange species and are not counted as synthesis targets.

Of the 101 target records, 100 were blocked in the zero-organic-pool startup
test. Examples with net certificates include ethanol, mannitol, limonene and
linalool. These are not new claims of physiological Cannabis production.
The certificates use 211 distinct balanced equations and preserve each equation's
full participants, coefficients, source records and candidate evidence IDs.
The longest selected solution uses 40 directed steps; this is not necessarily
the shortest pathway.

The numerical search uses SciPy/HiGHS linear programming
([solver documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html)).
It minimizes total nonnegative directed extent, reconstructs rational extents,
and rejects any solution that fails exact nondepletion or target-production
checks. Independent tests replay every published certificate with rational
arithmetic and verify carbon/charge conservation and CO₂-only net carbon input.
Solver-reported infeasibility remains numerical evidence, not an exact proof or
biological absence. Failed or incomplete numerical solves cannot become certificates.

**Pre-existing internal pools may be required.** Zero-net internal participants
are listed, but their origin, minimum pool quantities and feasible startup
sequence are not established. Closing a net balance does not trace atoms from
CO₂. The permissive carbon-free exchanges, unknown directions, possible
energy-generating cycles, absent thermodynamic constraints, and unresolved
compartments, specificity and activity prevent interpreting these results as
demonstrated pathways. Atom tracing remains deferred; the earlier startup and
finite execution certificates remain unchanged and separately labeled.

Install the optional numerical dependency with `pip install -e '.[flux]'` and
reproduce with `PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_net_flux`.
The report records the SciPy version and input checksums. The GCP snapshot is
`terpedia-489015.terpedia_core.cannabis_phase1_candidate_net_flux_20260904_v1`
(6,853 records: 6,220 targets, 100 certificates, 211 reactions, 321 compounds,
and one metadata record). The net-conversion explorer above displays this report
alongside, rather than replacing, the existing one-step map and evidence labels.

The [candidate-constrained scope and single-gap tests](data/phase1-candidate-scope.json)
separate chemistry-only connectivity from connectivity through equations with
attached candidate-enzyme evidence. Of 13,995 balanced equations, 1,472 have
candidate evidence and 12,523 are excluded in this diagnostic. An annotation or
homology candidate is not direct activity evidence. Both hypothetical reaction
directions remain allowed; no physiological reversibility claim is made.

- CO₂ alone: no additional compounds or targets become available.
- CO₂ plus all 101 carbon-free catalog species: nine additional compounds and
  one non-seed CannabisDB target (urea, CDB004840) become available. This leaves
  6,184 carbon-bearing target records blocked, one reachable and CO₂ itself seeded.

This is a **zero-organic-inventory startup diagnostic**, not a model of an
established plant's steady-state CO₂ assimilation. ATP, NAD(P)H, CoA and other
organic pools are not silently supplied. Failure to start from these seeds does
not demonstrate that a living plant cannot synthesize the metabolite. Conversely,
the permissive carbon-free reservoir and unknown directions mean even the urea
result does not establish an in-vivo Cannabis pathway. Atom tracing stays deferred.

There are 27 unsupported directed steps whose complete inputs are available and
whose outputs could expand this candidate-only scope. Each was tested separately
by temporarily admitting just its one full balanced equation (both hypothetical
directions), then repeating all-reactant expansion through candidate-linked
equations. Seven such tests add a target, never more than one per test. Across
these separate scenarios the additional target identities are uric acid,
carbon monoxide and hydrogen cyanide; they are not newly established pathways.
All 27 tests retain the full original source equation, prior screen outcome,
newly available compounds, target identities and complete prerequisite witnesses.
No hypothetical rescue is promoted into candidate-enzyme evidence. Curation must
first determine whether an equation is enzymatic, spontaneous, or a transformation
rule and assess its direction before deciding on protein-discovery work.

The reactant denominators are 1,626 unique reactants in the candidate-linked
subnetwork and 11,162 in the full network; neither uses inventory-only nodes.
All 12,440 target/scenario records, including blocked targets and seeds, remain
in the report. Reproduce with
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_candidate_scope`.
The GCP snapshot is
`terpedia-489015.terpedia_core.cannabis_phase1_candidate_scope_20260904_v1`
(12,719 records). Existing map evidence filters are unchanged by these
counterfactual tests; the report is a separate downloadable diagnostic.

The focused Cytoscape map now uses the
[combined enzyme overlay](data/phase1-combined-enzyme-overlay.json), joining
the earlier screen and the selected-route screen by exact balanced-equation ID.
The combined layer contains 744 screened-equation evidence records: 633 from
the earlier screen and 111 newly candidate-linked equations present in the
one-step map. The other 56 newly candidate-linked upstream equations are retained
in the [full route-screen overlay](data/phase1-route-enzyme-overlay.json), which
contains all 167 new evidence records. Every record links the correct full
screen report, all passing alignment IDs, representative alignments, reference
joins, candidate proteins and validation blockers. Representative alignments
are not functional rankings.

Candidate filters and edge details now cover 2,740 one-step hypotheses, 217 more
than the previous map; 17,032 hypotheses still lack candidate evidence. Target
coverage remains 289 carbon-bearing records with any candidate evidence: these
additions strengthen alternatives/upstream steps, not new target coverage.
All 19,772 one-step hypotheses retain their original blocked pathway status,
full coefficients, required inputs, target identities and direction assumptions.

The [route evidence status overlay](data/phase1-route-evidence-status.json) joins
the same evidence to all 304 complete prerequisite certificates without
rewriting their chemistry, finite quantities, seeds or execution assumptions.
All 304 routes still have missing candidate evidence; 343 distinct selected
equations remain without it. The first missing enzyme step is recomputed per
route. Neither this evidence update nor the preserved stoichiometric certificate
establishes a physiological Cannabis pathway. Atom tracing remains deferred.

Reproduce with `PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_route_overlay`,
then `PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_hypothesis_view`.
The GCP tables in `terpedia-489015.terpedia_core` are
`cannabis_phase1_route_enzyme_overlay_20260904_v1` (168 records),
`cannabis_phase1_combined_enzyme_overlay_20260904_v1` (745 records), and
`cannabis_phase1_route_evidence_status_20260904_v1` (305 records).

The [selected-route whole-proteome screen](data/phase1-route-protein-search.json)
compares all 30,304 pinned Cannabis proteins against all 11,487 retrieved
reference sequences (192 successful batches). DIAMOND sensitive search returned
34,100 alignments; 11,522 pass identity ≥30%, query coverage ≥50%, reference
coverage ≥50%, and E-value ≤10⁻⁵. These yield 619 distinct Cannabis candidates
and 1,197 protein–reaction hypotheses across 167 of the 510 selected-route gaps.
The remaining outcomes are 104 weak-hit-only equations, 118 no-hit equations,
and 121 without a reference sequence. Missing annotation or a failed homology
screen is not evidence that the reaction is biologically absent.

The newly candidate-linked equations occur in selected routes for 298 targets.
**All 304 routes still contain at least one equation without candidate evidence**
after this result is joined by exact balanced-equation ID. These counts do not
establish any complete Cannabis pathway. Direction, substrate specificity,
catalytic residues/domains, compartments and expression remain unverified.
Each result retains reaction-specific reference annotations, passing alignment
IDs, full reference and candidate sequences, exact search settings, source hashes,
weak-hit counts and proposed biochemical tests. The original route certificates
remain unchanged snapshots. The current map and route-status overlays above
integrate this screen without rewriting its source evidence or the certificates.

The verified GCP snapshot is
`terpedia-489015.terpedia_core.cannabis_phase1_route_protein_search_20260904_v1`
(24,331 records). Reproduce with:

```python
from pathlib import Path
from cannabis_carbon.phase1_new_protein_search import run
run(Path('data/reports/phase1-route-references.json'),
    Path('data/raw/phase1-route-protein-search'),
    Path('data/reports/phase1-route-protein-search.json'))
```

The [selected-route reference discovery](data/phase1-route-references.json)
retains all 510 enzyme-evidence gaps from the finite route certificates. Published
Rhea family mappings resolve 505 master families. Forty-three checksummed prior
lookup batches were reused and 16 new batches retrieved successfully. There are
11,487 distinct reviewed, nonfragment reference-protein leads spanning 389 of the
510 equations; the remaining 121 have no attached reference lead in this pass.
Exact source-ID matches are kept separate from direction-family matches. Reviewed
annotations are not necessarily direct assays and never establish Cannabis
activity or physiological direction. Full equations, route membership, lookup
URLs, snapshots and checksums remain attached. Reproduce with
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_route_references`.
The verified GCP snapshot is
`terpedia-489015.terpedia_core.cannabis_phase1_route_references_20260904_v1`
(12,057 records: 510 gaps, 11,487 references, 59 lookups and one metadata record).

The [route certificates and enzyme-gap queue](data/phase1-route-certificates.json)
expand all prerequisites for each of the 304 structurally reachable targets.
Each certificate retains full directed equations, exact rational reaction
extents, initial seed quantities, final inventory, and reaction-level candidate
evidence IDs. Independent replay requires every input quantity to be present
before any output is credited; every intermediate inventory stays nonnegative.
Tests additionally check total carbon and charge conservation. The only
carbon-containing initial reagent is exact CO₂. Atom tracing remains deferred.

These are deterministic first-witness routes, not shortest or optimal routes.
Conservative seed budgets ignore incidental coproduct credits when planning
upstream quantities; all coproducts are retained during replay. The 304 routes
use 740 distinct equations, of which 510 lack attached candidate-enzyme evidence.
Every selected route has at least one such gap. The largest certificate contains
172 reaction steps. The gap queue ranks equations by affected selected targets
and links their original reaction sources; it does not establish unavoidable
bottlenecks across alternative pathways. Candidate evidence is not direct activity.

Finite sequential stoichiometric feasibility is stronger than graph adjacency,
but it does not validate physiological direction, thermodynamics, compartments,
enzyme activity, or the highly permissive carbon-free seed reservoir. Every
route remains biologically unestablished. All 12,440 target/scenario records are
preserved, including blocked targets and explicit seeds. The full network and
scope reports are pinned by checksum. Reproduce with
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_routes`.
The report and export are published to
`terpedia-489015.terpedia_core.cannabis_phase1_route_certificates_20260904_v1`:
12,440 target records, 304 routes, 510 enzyme gaps, and one metadata record.
The static map remains a one-step view; this release publishes downloadable
multi-step certificates rather than a route-explorer UI.

The [all-reactant scope audit](data/phase1-all-reactants-scope.json) moves beyond
single-edge participation: a reaction can expand the available-compound set only
when **every** input is already available. It uses the
[full balanced network](data/phase1-full-balanced-network.json), which adds 3,120
upstream equations outside the earlier target-matching subset. The resulting
network retains 13,995 balanced full equations, 16,938 exact compound structures,
all 6,220 target records, source orientations, and candidate-evidence links.
All 36,552 Rhea source records were considered: 28,102 pass balance, 8,408 are not
auditable with concrete structures, and 42 are imbalanced. Reverse encodings are
deduplicated without canceling or rescaling participants.

Two **explicitly non-physiological structural scenarios** are retained:

- **CO₂ alone:** one seed, no newly available compounds; 6,185 carbon-bearing
  target records remain blocked and the CO₂ target itself is a seed.
- **CO₂ plus every cataloged carbon-free species:** 102 seeds (CO₂ plus 101
  carbon-free species), 6,486 newly available compounds, and 304 structurally
  reachable carbon-bearing targets; 5,881 carbon-bearing targets remain blocked.
  This permissive reservoir includes unusual inorganic/redox species and is **not
  a proposed plant growth medium**. No other carbon-containing compound is seeded.

Every reaction is tested in both directions only as an upper-bound diagnostic,
not an assertion of physiological reversibility. All input coefficients remain
in the full equations and witness records. Qualitative availability does not
model finite quantities, flux bounds, thermodynamics, compartment compatibility,
or enzyme activity. Each non-seed witness uses inputs from strictly earlier
expansion levels; missing cofactors and unseeded cycles cannot bootstrap themselves
through graph adjacency. The deepest witness reaches 97 expansion levels, not
necessarily 97 distinct reactions. The denominator for reactant coverage is the
11,162 unique reactants in the all-directions scenario, not all inventory nodes.
This result is **not** a claim that 304 Cannabis pathways are established, and
atom tracing remains deferred.

GCP snapshots:
`terpedia-489015.terpedia_core.cannabis_phase1_full_balanced_network_20260904_v1`
(45,604 records) and
`terpedia-489015.terpedia_core.cannabis_phase1_all_reactants_scope_20260904_v1`
(19,032 records). Full source/checksum metadata, both seed lists, all targets,
missing inputs, and step witnesses are retained. Reproduce using
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_scope`.
The focused map remains a one-step view; these multi-step witnesses are available
in the linked report pending dedicated route-view integration.

The focused Cytoscape map now includes the
[screened-enzyme evidence layer](data/phase1-screened-enzyme-overlay.json).
It joins 633 independently balanced equation IDs to the new Cannabis homology
screen without altering source equations, coefficients, chemical identities,
directions or pathway status. Candidate evidence now covers **289 carbon-bearing
targets**, up from 79; 210 additional targets gain candidate leads. Across all
target projections, 1,487 hypotheses gain attached evidence, bringing the total
to 2,523 (including carbon-free targets). These are not confirmed enzyme counts.

The map's candidate/missing-evidence filter uses this combined snapshot.
Each added evidence record retains all passing alignment IDs and an explicitly
labeled representative alignment per protein, with identity, both coverage values,
reference accession, and unresolved-specificity/direction warnings.
The base hypothesis catalog remains immutable; the evidence layer has its own
parent/search checksums and GCP snapshot:
`terpedia-489015.terpedia_core.cannabis_phase1_screened_enzyme_overlay_20260904_v1`.
Its 634 rows contain 633 evidence records plus integration metadata. Rebuild using
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_screened_overlay`, then
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_hypothesis_view`.
All 19,772 hypotheses remain blocked pending biochemical and pathway validation;
the 409 carbon-bearing targets with structural production hypotheses are unchanged.

The [new reference-driven Cannabis proteome screen](data/phase1-new-protein-search.json)
retrieved and validated all 4,653 reference sequences, then searched all 30,304
proteins in the checksum-verified Cannabis reference proteome using DIAMOND
sensitive mode with unlimited reported targets and one HSP per pair. Of 56,294
pairwise alignments, 17,407 pass the declared screen: identity at least 30%, both
query and reference coverage at least 50%, and E-value at most 1e-5.

The screen identifies **1,074 distinct Cannabis candidate proteins**, producing
8,651 protein–reaction hypotheses across 633 balanced equation gaps. Of the 330
priority carbon-bearing targets previously lacking any candidate-supported
alternative, **210** now have screened homology candidates. This is not experimental
enzyme confirmation or increased CO₂ reachability. Reaction direction, exact
substrate specificity, catalytic residues/domains, tissue/compartment compatibility,
and all-input availability remain unresolved. All 2,567 original equation gaps
are retained, including missing-reference, weak-only, and no-hit outcomes.

The report embeds passing alignments once, with per-equation alignment IDs,
reference annotations, complete candidate/reference sequences and hashes, retrieval
URLs, the exact search command/version, and validation blockers. The complete raw
alignment output (including weak hits) is checksummed separately; per-equation raw
counts distinguish weak evidence from no hits. The focused map incorporates these
candidates through the separate, checksummed screened-enzyme layer described above.

GCP snapshot:
`terpedia-489015.terpedia_core.cannabis_phase1_new_protein_search_20260904_v1`.
Its 25,780 records include equation gaps, reference sequences, Cannabis candidates,
passing alignments, retrieval provenance, and search metadata. Reproduce via
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_new_protein_search`.
Atom tracing remains deferred.

The [new carbon-target enzyme-reference discovery](data/phase1-new-references.json)
retains all 3,000 carbon-bearing one-step hypotheses lacking attached candidate
enzyme evidence, grouped into 2,567 balanced equations. It prioritizes the 330
targets with no candidate-supported alternative. Queries covered 1,188 explicit
Rhea master families in 48 successful UniProt batches, restricted to reviewed,
non-fragment records. Returned annotations are joined through the checksum-verified
published Rhea direction-family table, never guessed from adjacent numeric IDs.

Reference leads were found for 1,140 gap equations and 276 of the 330 priority
targets, involving 4,653 distinct reference proteins. These are **reference
annotations, not newly identified Cannabis proteins**. Exact source-ID matches
and family-only matches remain separate; neither establishes the producing
direction required by a hypothesis. Unqueried alternative gaps and queries
without matching references are retained explicitly. Next: validate reference
sequences, screen the full Cannabis proteome, and evaluate substrate specificity.
The live graph's candidate-enzyme counts are not inflated by these references.

GCP snapshot:
`terpedia-489015.terpedia_core.cannabis_phase1_new_references_20260904_v1`.
Reaction-gap rows, reference protein annotations, lookup URLs and snapshot
checksums are retained in `data/reports/phase1-new-references.ndjson` and the
published report. Regenerate via
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_new_references`.

The [interactive balanced-hypothesis map](hypotheses.html) now exposes this
catalog directly in Cytoscape. Search by label or CannabisDB ID, filter by
attached candidate enzyme evidence, and select one hypothesis. All 6,220 records
remain searchable, including explicit gaps. Each selected graph includes every
input and output, with compounds as nodes and directed reaction projections as
edges. Multiple projected edges represent **one** full reaction, not separate
pathway steps or traced carbon flows. The equation, coefficients, source records,
enzyme evidence, blockers, bootstrap warnings, and proposed tests remain inspectable.

The view is static GitHub Pages content, derived from the same GCP-backed catalog.
It loads a small target list and one reaction shard on demand; every reaction
shard is below 0.6 MB. `data/hypothesis-view/index.json` preserves the parent
report checksum and all generated-file hashes. Regenerate the view using
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_hypothesis_view`.
The whole-network map remains available and links to this focused view.

The [full-inventory one-step hypothesis catalog](data/phase1-target-hypotheses.json)
distinguishes reaction participation from **net production**. Of 6,220 CannabisDB
records, 424 have balanced net-production hypotheses, six occur only unchanged
on both sides of balanced equations, six match only balance-unresolved equations,
and 5,784 have no exact structure match. The unchanged group is aldehydo-L-rhamnose
plus five inorganic elements; their occurrence does not establish a producing step.
Transport and compartment distinctions are not inferred from identical structures.

The carbon-specific denominator is **6,186 carbon-bearing records**, with **409**
having net-production hypotheses; 34 carbon-free records remain explicitly tracked
for nutrient uptake, transport, and input supply. There are 3,231 carbon-bearing
target hypotheses out of 19,772 hypotheses overall. Water and ion reactions dominate
the latter count, which must not be presented as carbon-pathway completeness.

The catalog retains 15,669 exact compound structures, including every resolved
CannabisDB target and every participant in the retained balanced equations.
Its 10,875 full-equation records deduplicate exact equations and reversed encodings
from 22,828 source records while retaining each source orientation and provenance.
Coefficients and participants are not canceled away or rescaled for deduplication.
Every hypothesis includes full inputs/outputs, a positive net target coefficient,
source-linked candidate enzyme evidence, explicit blockers, and proposed tests.
There are 1,036 hypotheses with attached candidate enzyme evidence and 18,736
without it; neither count establishes enzyme activity. Nine hypotheses require
the target itself as an input and are flagged as needing a bootstrap supply.
All hypotheses remain blocked: physiological direction, all-input supply,
compartment compatibility, and Cannabis enzyme activity are unconfirmed.
Both directions are examined only as a structural hypothesis scenario. No seeds,
including intracellular currencies, are silently provided; atom tracing is deferred.

GCP snapshot:
`terpedia-489015.terpedia_core.cannabis_phase1_target_hypotheses_20260904_v1`.
Regenerate with
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_target_hypotheses`;
the module's `export_table` function produces the local NDJSON load artifact.
The roughly 65 MiB static download is separate from the Cytoscape working map.

The [full-Terpedia-Rhea participation audit](data/phase1-target-rhea-coverage.json)
extends the working-network audit below with all 36,552 Rhea source equations
from the checksum-verified Terpedia catalog snapshot. It retains all 6,220
CannabisDB records: **430** match participants in balanced equations (**321 new**
relative to the working network), six match only balance-unresolved equations,
and 5,784 have no exact encoded-structure match. New leads include eugenol,
nerol, methyl salicylate, xylitol, and vanillic acid. These are source-database
reaction leads, not evidence of Cannabis enzyme activity or complete CO₂ routes.

Across the overlapping core, expansion, and full-Rhea source layers, 38,150
scoped equation records were independently audited: 29,310 balanced, 42
imbalanced, and 8,798 not auditable with concrete structures. These are **not
deduplicated biochemical reaction counts**. The downloadable report retains
the 29,538 equations with at least one exact target participant, preserving full
source equations and all required participants. A right-side match is not a
physiological producer assignment; all-input availability remains unestablished.
Unspecified stereochemistry is flagged and never upgraded to an exact known
stereoisomer. Atom tracing remains deferred.

GCP snapshot:
`terpedia-489015.terpedia_core.cannabis_phase1_target_rhea_coverage_20260904_v1`.
It contains 6,220 target rows and 29,538 matching equation rows. Regenerate the
report and local NDJSON load file using
`PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_target_rhea_coverage`.
The large derived NDJSON is ignored by Git; source checksums and the GCP load
receipt are retained. The static report is approximately 58 MiB and is not
automatically loaded into Cytoscape or silently promoted into its working map.

The [whole-CannabisDB strict participation audit](data/phase1-target-coverage.json)
retains all 6,220 source records (6,203 distinct canonical encoded structures).
Against the core NetworkDB and balanced expansion catalog, 109 records match
participants in independently balanced equations, four match only balance-unresolved
equations, and 6,107 have no exact encoded-structure match in these reaction sets.
Two of the 109 matched records have unspecified/unknown potential stereochemistry,
which remains flagged. Matching does not relax protonation, salts, isotopes,
tautomers or stereochemical encoding. Twelve additional records have only
source-supplied ChEBI reaction links; those are not structural confirmation.
These are participation and data-linking metrics, not demonstrated biosynthetic
coverage or evidence that unmatched compounds cannot be made by Cannabis.

The audit preserves 1,598 scoped reaction records, not necessarily distinct
biochemical reactions: 717 core equations pass a concrete structure-based balance
audit, 390 core equations are not auditable by that method, and all 491 expansion
catalog equations pass. Generic participants are not labeled imbalanced merely
because their concrete formulas are unavailable. Equation left/right roles are
not automatically physiological direction. No CO₂ reachability is asserted.

GCP table:
`terpedia-489015.terpedia_core.cannabis_phase1_target_coverage_20260904_v1`.
Its 7,818 records retain both the 6,220 targets and 1,598 source reaction audits,
with source checksums, structure flags, roles and coefficients, xref-only leads,
and explicit unmatched outcomes. Load source:
`data/reports/phase1-target-coverage.ndjson`.

The [consolidated Phase 1 reaction catalog](data/phase1-reaction-catalog.json)
contains 491 balanced source-reaction variants, 215 exact compound structures,
and an audit ledger retaining all 690 original expansion variants. Reactions
store full reactant/product lists with coefficients and include carbon-free
co-reactants. Structure IDs preserve canonical isomeric SMILES identity; no
CannabisDB cross-reference is inferred from those IDs. Source-linked homology,
core annotations and direction-unresolved family evidence remain separate.

The 108 balanced-alternative links resolve to only 37 distinct Rhea equations,
all already among the 491 balanced variants. They are not 108 new reactions.
The catalog attaches each link and participant difference without automatically
transferring the original MARTS protein annotation to the Rhea equation. All 199
non-balanced originals remain excluded from the balanced set; 108 have reference
alternatives and 91 do not. This catalog covers the expansion, not the complete
Cannabis metabolome or core NetworkDB, and does not establish CO₂ reachability.

GCP catalog table:
`terpedia-489015.terpedia_core.cannabis_phase1_reaction_catalog_20260904_v1`.
Its 1,396 rows distinguish `compound`, `reaction` and `source_variant` records,
with stable IDs, complete record JSON and a report checksum. Load source:
`data/reports/phase1-reaction-catalog.ndjson`.

The [source-backed balance alternatives](data/phase1-balance-reference.json)
compare all 199 non-balanced expansion variants against 36,552 Rhea equations
in Terpedia's normalized reaction catalog. Matching preserves stereochemistry,
isotopes, formal charge and the carbon-containing species on each reaction side.
There are 108 links from imbalanced MARTS variants to 37 distinct balanced equations:
96 preserve carbon-compound coefficients, while 12 require coefficient changes.
Every alternative passes an element-and-charge audit and explicitly records
added/removed participants and coefficient differences. Original records remain
unchanged. These source equations are candidate alternatives, not automatic
confirmation of MARTS enzyme associations or physiological direction.
The remaining 42 imbalanced and 49 generic/unauditable variants have no exact
concrete source alternative in this audit.

GCP snapshot:
`terpedia-489015.terpedia_core.cannabis_phase1_balance_reference_20260904_v1`.
It retains all 199 gaps with original equations, reference equations, structured
participants, balance results, differences and provenance; local load source:
`data/reports/phase1-balance-reference.ndjson`.

The separate [direction-unresolved family search](data/phase1-family-protein-search.json)
searched all 30,304 Cannabis proteome sequences against all 168 retrieved reviewed
reference proteins using DIAMOND's sensitive mode. Among 282 retained gap variants,
62 have screened family-homology candidates (172 distinct Cannabis proteins;
831 protein–reaction hypotheses), 41 have weak hits only, 22 have no hits, and
157 have no family-reference sequence to search. Screening requires at least
30% identity, 50% coverage of both proteins, and E-value at most 1e-5.
Full alignments, candidate/reference sequences and hashes, the search command,
source reaction-family annotations, and proposed validation tests are retained.
All requested reaction directions remain unverified; this layer is not silently
added to the main map's 202/491 candidate-evidence metric.

The matching GCP snapshot is
`terpedia-489015.terpedia_core.cannabis_phase1_family_protein_search_20260904_v1`,
with one row per gap variant and JSON columns preserving alignments, references,
candidate proteins and validation blockers. The load source is
`data/reports/phase1-family-protein-search.ndjson`.

The [Phase 1 reference-discovery audit](data/phase1-reference-discovery.json)
checks all 282 missing-reference variants against Terpedia's current normalized
reaction catalog (282 exact reaction-ID/SMARTS matches) and reviewed, nonfragment
UniProt entries. No returned annotation matches the exact requested directional
Rhea ID. Explicit joins through the published Rhea direction-family table provide
168 reference proteins for 125 variants. These are direction-unresolved search
inputs, not newly screened Cannabis candidates, so the map's candidate-enzyme
completeness counts remain unchanged.

GCP snapshot table:
`terpedia-489015.terpedia_core.cannabis_phase1_reference_discovery_20260904_v1`.
It contains one row per gap variant, with separate exact-reference and
direction-family-reference JSON columns, source catalog rows, direction mappings,
sequence-search status and checksums. The reproducible local load source is
`data/reports/phase1-reference-discovery.ndjson`; the full audit retains query URLs
and raw snapshot checksums. The discovery module uses `bq` on PATH, or the
`CANNABIS_BQ` environment override.

Reference interpretation follows the [Rhea API documentation](https://www.rhea-db.org/help/rest-api)
and its [published directional mapping](https://www.rhea-db.org/help/download).
Reviewed annotations and direction-family membership do not establish a
direction-specific Cannabis reaction or a complete pathway.

The static map's Phase 1 evidence filter uses
`docs/data/phase1-map-evidence.json`, joining exact reaction-ID/SMARTS variants
without a name or reaction-family fallback. Its expansion-only denominator is
690 variants: 491 balanced, 150 imbalanced, and 49 not auditable. Of the balanced
variants, 178 have screened homology candidates and 24 additional variants have
core enzyme associations; 289 lack either evidence type (282 missing references,
four no-hit searches, three weak-hit-only searches). These are not confirmed
enzyme or complete-pathway percentages. The filter dims nonmatching edges and
nodes, composes with existing filters, and updates when the candidate layer loads.

The Phase 1 core enzyme evidence supplement
(`docs/data/phase1-core-enzyme-evidence.json`) links 24 balanced expansion
variants lacking search reference sequences to existing core enzyme associations
through explicit Rhea IDs or listed directional IDs. It preserves association
sources, qualifiers, core equation participants, and input checksums. These are
annotation/candidate links, not new biochemical confirmations; original sequence
search statuses remain unchanged, and core equation orientation does not establish
the expansion's physiological direction.

The Phase 1 experimental shortlist (`docs/data/phase1-experimental-shortlist.json`)
groups 3,128 protein–reaction proposals under 18 Cannabis proteins, covering
178 balanced reaction variants. Every protein retains its source FASTA header,
sequence and checksum; each proposal includes ranked reference alignments,
reaction SMARTS, provenance and unresolved validation steps. The proposals
share homology evidence and must not be counted as independent confirmations.

Targeted Phase 1 protein search: all 123 source references were resolved to
sequences, including 97 UniProtKB, eight UniParc, 16 NCBI protein records and
two single-CDS translations annotated in NCBI nucleotide records. The latter
preserve separate identities: MK803261.1 → QDZ36304.1 and
MK803262.1 → QDZ36305.1, with source XML checksums and CDS locations.
Ambiguous multi-CDS records are rejected, not assigned an arbitrary protein.
DIAMOND searched all 30,304 proteins in UP000583929. Of 491 balanced reaction
variants, 181 have hits and 178 have candidates passing 30% identity and 50%
coverage of both query and reference (E-value at most 1e-5). These correspond
to 18 distinct Cannabis proteins. The report at
`docs/data/phase1-targeted-protein-search.json` retains every balanced variant,
no-hit outcomes, missing references, sequence hashes and search command.
Homology supports experimental prioritization, not established activity.

NCBI reference WP_169336908.1 (source-annotated squalene/oxidosqualene cyclase
from Eudoraea adriatica) supplied the previously unsearched reference for
MARTS:6816, MARTS:6817 and MARTS:6818. Each now has ten screened Cannabis
protein candidates, adding 30 reaction–protein hypotheses, not confirmed enzymes.

Phase 1 priority: atom tracing is deferred at the user's request. The
[enzyme discovery queue](data/phase1-enzyme-discovery-queue.json) groups 690
reaction-ID/SMARTS variants: 185 balanced variants have source protein IDs,
306 balanced variants need characterized reference enzymes, and 199 require
stoichiometry review. Source protein IDs do not establish Cannabis activity.

Atom-continuity correction: traversal now uses compound and atom index, finding
continuous CO₂-to-core witnesses for 2,308 candidate bridge records; 1,326
retain unresolved entity-path fallbacks. All 3,634 still have unresolved
end-to-end carbon provenance pending core-to-target identity mapping and
required-input verification. Reverse traversal swaps both entity
and atom endpoints. The 2,154 complete product-carbon mappings describe local
identity-pair correspondences only; they do not establish continuous CO₂
provenance. Balance eligibility checks only the candidate bridge reaction.

Generated by `cannabis_carbon completeness`:

- Terpedia metabolites without reaction participation: **7**.
- Terpedia has **1,107** working reactions, including **9** explicitly
  non-enzymatic reactions and **1,098** enzyme-requiring reactions.
- **12** enzyme-requiring reactions lack a direct enzyme association; all **12**
  have attached candidate proteins.
- The completeness artifact now separately reports the expanded hypothesis
  layer: **7,384** directed edges (**6,429 candidate**, **955 unresolved**),
  **5,429** source reaction hypotheses, **1,324** product records, and **1,121**
  edges with either enzyme-gene or curated enzyme-catalog evidence. These
  counts do not alter the core balance or CO₂-lineage denominators.
- CannabisDB compounds: **6,220**.
- CannabisDB carbon atoms: **285,623**.
- PubChem resolution: **3,586 exact**, **130 ambiguous**, **10 connectivity
  candidates**, **180 name candidates**, and **2,314 explicit no-matches**; all
  6,220 records have a valid InChIKey. Connectivity and name matches remain
  candidate identity evidence.
- CannabisDB XML external-ID coverage: **3,737 records** have at least one
  first-party external identifier; **2,483** have none. These source links are
  retained separately from exact RDKit identity matches.
- Terpedia reaction-product carbon atoms: **22,751**.
- Structurally inferred reaction-product carbon mappings: **5,225** (**22.966%**).
- Additional reaction-product carbons are explicitly classified as **8,720
  candidate**, **8,440 ambiguous**, and **366 unresolved**; every product
  carbon has a mapping row.
- RDKit MCS mapping uses one documented **10-second per-pair timeout** so
  repeated report generation is reproducible; timeout failures remain
  explicitly unresolved rather than silently assigned.
- CannabisDB carbon coverage remains **not computable** until a verified
  complete pathway crosswalk is available.
- The consolidated candidate CO₂-path carbon layer contains **3,634**
  reversible-upper-bound paths covering **79** candidate products. **2,154**
  paths have complete product-carbon mappings, with **132,378** mapped and
  **4,220** unresolved product carbon atoms; **22,649** atom-level core-path
  edges are retained. These are structurally mapped hypotheses, not proof of
  physiological direction or in-vivo biosynthesis.
- The candidate path rows now carry the source SMARTS balance join: **1,692**
  paths are balance-eligible, **104** are imbalanced, and **1,838** are not
  auditable. Non-eligible rows remain retained for review and are not promoted
  into Phase 1-valid pathway steps.
- Candidate-path target lineage now retains explicit CannabisDB IDs for **7**
  candidate compounds spanning **230** CannabisDB carbon atoms; these are
  connectivity-candidate labels, not exact identity or biosynthesis claims.
- The candidate expansion balance audit covers **32,904** source edges and
  **690** unique reaction/SMARTS records: **30,045** edges are balanced,
  **760** imbalanced, and **2,099** not auditable. At the unique-reaction
  level, **491** are balanced, **150** imbalanced, and **49** not auditable.
  The report is available as `data/reports/terpene-identity-set-candidate-
  expansion-balance-audit.json` and the published copy under `docs/data/`.
- CannabisDB-to-Terpedia identity crosswalk: **152 exact**, **1 ambiguous**,
  **6,068 without an exact identity**, and **5,769 without any identity
  resolution**; 1,585 CannabisDB carbon atoms are identity-linked. The
  reverse-direction Terpedia inventory has 1,113 metabolites without a
  CannabisDB match. A separate connectivity candidate layer contains 250
  one-to-one candidates and 40 ambiguous groups; the candidate-only RDKit
  canonical-tautomer layer currently adds 0 one-to-one candidates and 1
  ambiguous group. The unique-name candidate layer currently adds 8
  one-to-one candidates and is retained as weak, non-exact identity evidence.
  The CannabisDB XML ChEBI identifier layer adds **100 structure-verified
  identifier links** and records **6 identifier conflicts** separately; these
  links are not merged into the exact InChIKey count. These
  alternatives may differ in stereochemistry, protonation, or tautomer state
  and are never treated as exact identity.
- The BigQuery table `terpedia-489015.terpedia_core.terpene_identity_set` was
  queried by full InChIKey and by the first (connectivity) InChIKey block. The
  live table currently contains 268,924 identity rows (199,234 with an
  InChIKey); the checked-in retrieval contributes **246 exact identity-set
  matches** plus **104 connectivity-only candidate compounds** across **437
  candidate identity records**
  to CannabisDB, with source memberships, source record IDs, manifest URIs,
  and source-file URIs preserved in
  `data/reports/terpene-identity-set-match.json`.
- The separate
  `data/reports/terpene-identity-set-connectivity-upstream.json` queue joins
  connectivity-only identity candidates to **2,463** directly characterized
  MARTS-DB upstream edges covering **74** candidate products and **31**
  reactions. Candidate CannabisDB IDs and all identity alternatives are
  retained; these records are producer hypotheses, not exact identities or
  confirmed Cannabis pathways.
- A separate snapshot of
  `terpedia-489015.terpedia_core.terpene_biotransformation_hypothetical_connections_current_v2`
  preserves **193** current source-linked candidate edges in
  `data/terpedia/terpene-biohypothetical-connections-current-v2.json`. After
  normalizing legacy `T000`/`T100` identity aliases, all **193** rows match
  edges already present in the 1,901-edge hypothesis layer; the current table
  adds no new edge, but provides an independently versioned GCP cross-check.
- The exact identity-set match now has a dedicated producer-gap queue at
  `data/reports/identity-set-pathway-gap-queue.json`: **161** of 246 exact
  identity-set targets do not appear as products in the refreshed upstream
  edge extraction, representing **2,551 carbon atoms**. Ten already participate
  in a core reaction record; all remain queued for producer-reaction and
  Cannabis-proteome searches because identity alone is not pathway evidence.
- The live curated BRENDA/Rhea SMARTS table contains **86** structure-resolved
  reaction templates (77 terpene-biosynthesis or rearrangement roles). The
  snapshot is preserved in
  `data/terpedia/terpene-brenda-reaction-smarts-curated-current-v2.json`.
  RDKit found no direct product-structure matches among the 6,220 CannabisDB
  compounds, so these templates remain enzyme/reaction evidence and do not
  create unsupported CannabisDB pathway edges.
- CannabisDB-wide pathway coverage remains **not computable** because identity
  linkage is incomplete and mapped reaction products are not yet connected to
  CO₂ by a complete directed reachability proof.

The import now preserves CannabisDB source names and aliases in
`docs/data/compounds.json`. The generated
`data/reports/named-specialty-inventory.json` identifies 88 records (1,910
carbon atoms) by an explicit cannabinoid/cannabis specialty-name selection
rule. Sixteen have exact Terpedia identities and explicit reaction
participation, including the CBGA oxidocyclization and acid-decarboxylation
edges. Its v2 `review_queue` ranks all 88 records by unresolved carbon atoms
and retains PubChem status, the lineage blocker, and the reversible structural
upper bound for each record. The remaining records are documented curation
gaps, not an inference that the compounds are absent from Cannabis or that
they share one route.

## CO₂-only carbon lineage

The report `data/reports/carbon-lineage.json` applies the plant carbon-source
constraint explicitly: CO₂ is the only admissible carbon input. In the current
working network, 2,641 carbon atom nodes across the reachable metabolite
entities are reachable from the CO₂ seed through 3,990 inferred structural
edges and 54,509 explicitly retained candidate alternative edges. Candidate alternatives
come from ambiguous one-to-one RDKit mappings and pairwise MCS candidates for
small multi-substrate reactions; they are visible for review but are not
promoted to confirmed provenance. This uses 22 explicit Terpedia
physiological direction assertions plus source-backed Rubisco and Calvin-cycle
direction records. There are 460 carbon-containing reactant entities outside
that lineage and therefore reported as external-carbon-source blockers. The
current target summary is 1 supported, 199 candidate, and 6,020 unresolved
CannabisDB compounds. The one supported target is the CannabisDB CO₂ record
itself; 199 candidate targets have exact or candidate identities with partial
or complete CO₂-reachable carbon atoms. The remaining connectivity identity
candidates remain unresolved because they are not yet reachable from CO₂. This
indicates that the imported reaction network is
missing CO₂-assimilation and/or structure-resolved links needed to establish
the plant-wide carbon path; it is not evidence that Cannabis uses another
carbon source.

The completeness artifact now includes `co2_lineage.target_triage`: 5,795
targets have no Terpedia identity, 225 have an exact or candidate identity but
remain unresolved for CO₂ reachability, and 200 identity-resolved targets are
already candidate or supported. Of the 88 specialty-name targets, 17 are
candidate and 71 remain unresolved. This partitions the next work into identity
resolution versus reaction/pathway reconstruction.
It also reports `carbon_mapping_blockers`: 493 reactions and 8,806 product
carbon rows currently require ambiguity or unresolved-mapping curation.

The companion `data/reports/carbon-atom-audit.json` and published
`docs/data/carbon-atom-audit.json` artifacts partition all **285,623**
CannabisDB carbon atoms into **1 supported**, **1,421 candidate**, and
**284,201 unresolved** atoms. Each group retains CannabisDB atom indices,
resolved Terpedia atom indices when available, reaction references, provenance
URLs, and the blocking reason. The target atom indices are explicitly scoped to
RDKit atom indices in each CannabisDB SMILES field and are not assumed to equal
the source-SDF atom ordering. The groups are required to partition the complete
carbon-atom set for each compound. For every reachable atom, `co2_paths` now
serializes one CO₂-to-atom chain with entity-local atom indices, reaction IDs,
edge evidence status, provenance, known enzyme IDs, and genome candidate-protein
IDs from NetworkDB; unresolved atoms have no fabricated path and retain their
blocker. These enzyme and protein annotations remain hypotheses or source
associations, not functional validation.
The completeness artifact reports the same denominator as percentages:
0.000350% strict-supported, 0.475102% evidence-bearing when candidate atoms
are included, and 99.524898% unresolved.

These are data-coverage metrics. They do not imply that all CannabisDB
compounds are endogenous cannabis metabolites or that Terpedia's reaction graph
constitutes a demonstrated in-vivo pathway.

## Unified NetworkDB

The Pages visualization loads the complete `docs/data/network-map.json` by
default: 10,148 compounds, 1,107 reaction records, 7,384 general hypothesis
edges, and 2,495 distinct connectivity-candidate producer edges. The latter
are rendered as a separate purple dashed layer and remain candidate-only.
The compact CO₂-focused projection remains available at
`docs/data/network-map-focus.json` (1,261 compounds, 65 reaction records, and
3,450 general hypothesis edges) and through `?full=0`; the browser fallback graph is
not the published map. Selecting a compound in the map now reports how many of its carbon atoms
have explicit serialized CO₂ paths and the maximum number of reaction steps;
the corresponding atom-by-atom chains remain in the downloadable audit.
The separate `docs/data/hypothesis-lineage.json` artifact traverses only
candidate hypothesis edges from the core CO₂-reachable set. It currently adds
5 provisional CannabisDB targets (50 carbon atoms); all 5 targets
(50 carbon atoms) have complete candidate atom mappings across their
hypothesis paths. The remaining 955 unresolved edges are excluded because
their substrate structures are missing.

The complete CannabisDB XML export is preserved at
`data/terpedia/cannabisdb-compounds.xml.gz` with its retrieval URL and SHA-256
in `data/terpedia/cannabisdb-source-manifest.json`. The normalized table at
`data/terpedia/cannabisdb-compounds.json` contains all 6,220 accessions and
retains source structures, descriptions, synonyms, external identifiers, and
general references. The XML export contains first-party cross-references to
PubChem, ChEBI, KEGG, ChemSpider, FoodDB, BioCyc, and other databases; these
are now extracted rather than discarded by the SDF-only importer.

`docs/data/networkdb.json` is the published normalized inventory used for
cross-source reconciliation. It contains 6,220 CannabisDB compound records,
1,267 working Terpedia metabolite records, 1,107 working biochemical reaction records,
all reaction participants and stoichiometric coefficients, enzyme associations,
source URLs, and the 152 exact CannabisDB–Terpedia identity links. The two
namespaces remain separate so an identity match cannot silently collapse distinct
source records; links are explicit in `identity_links` and `identity_link`.
NetworkDB also carries 273 connectivity-level candidate identity links and 41
ambiguous candidate groups as reviewable, non-exact alternatives, including
the canonical-tautomer candidate layer. These links do not override exact
identity.
Records with a matching Terpedia identity-set row carry that row's
source-linked identity evidence in `terpedia_identity_set`; the current
snapshot contains 246 exact matches. Connectivity-only alternatives are
retained separately in `terpedia_identity_set_candidates` and are never
merged into exact identity. This table is an identity resource, not evidence
of endogenous Cannabis biosynthesis or pathway direction.
The NetworkDB field `identity_candidate_upstream_connections` adds 2,495
source-linked producer hypotheses generated from the 2,463-row MARTS-DB queue;
multiple CannabisDB candidates are retained when one connectivity prefix maps
to more than one compound. These edges are intentionally excluded from the
balanced reaction and CO₂-lineage counts.
The companion `data/terpene-identity-set-candidate-expansion.json` expands the
37 distinct candidate precursors through three upstream Terpedia hops, yielding
32,904 source-linked edges across 729 products and 1,032 precursors (18,780
Rhea and 14,124 MARTS-DB records). It is a deeper hypothesis queue, not part of
the confirmed reaction or CO₂-lineage counts.
The bridge audit `data/terpene-identity-set-candidate-expansion-bridges.json`
finds 15,266 structure-preserving links from this neighborhood to 117 product
and 147 precursor identities matching Terpedia core structures. None currently
touches the directed CO₂-reachable core, so these branches remain explicitly
blocked at carbon-source connection. Under the separate all-reactions-reversible
sensitivity run, 2,308 bridge records touch the reversible CO₂ upper bound across
90 reaction records; this identifies reaction-direction curation as a concrete
next discriminator without claiming physiological reversibility.
The companion `data/terpene-identity-set-candidate-expansion-carbon-mapping.json`
applies RDKit identity-pair carbon correspondence to 6,607 unique bridge pairs:
251,320 product-carbon correspondences are mapped, 5,143 pairs are fully
inferred, and 4,060 product-carbon rows remain unresolved. These mappings are
candidate structural evidence and are not added to the directed CO₂ lineage.
The reversible sensitivity artifact
`data/terpene-identity-set-reversible-candidate-lineage.json` enumerates 3,634
ordered candidate paths from the CO₂ seed to 79 candidate products through 7
core anchors. These paths are useful for direction and enzyme experiments but
are explicitly not directed biological pathways.
The merged carbon-provenance artifact
`data/terpene-identity-set-reversible-candidate-lineage-carbon.json` joins
those paths to the RDKit pair mappings: 2,154 paths have complete product-
carbon correspondence, covering 132,378 mapped product carbons, while 4,220
product-carbon rows remain unresolved. These are reversible sensitivity
hypotheses now retain atom-level CO₂-to-core carbon edges and do not alter the
directed atom audit.
When available, each CannabisDB record also carries the exact-InChIKey PubChem
resolution status and CID/structure properties; unresolved PubChem queries are
retained as explicit negatives rather than being treated as missing data.
The focused upstream evidence queue
(`data/terpedia-identity-set-upstream-mapped.json`) joins current Terpedia
metabolic-map edges to exact identity-set matches and runs the RDKit carbon
mapper on each reaction SMARTS. It currently contains 3,401 edges, 87
identity-set target structures, 1,373 precursor identities, and 364 reactions;
1,530 edges are candidate-mapped, 874 are inferred, and 997 remain unresolved,
with 29,547 product-carbon atoms explicitly unresolved in whole-reaction
mapping. It remains separate from
the balanced reaction catalog and CO₂ lineage until each edge receives
stoichiometric, direction, enzyme, and mapping review.
The same report carries a separate target-pair mapping: 2,567 of the 3,401
precursor/product identity pairs have complete inferred or candidate
carbon-skeleton mappings, while 834 remain unresolved; 77,640 target-pair
carbon atoms are mapped and 2,297 remain unresolved. This pair layer is not substituted for
whole-reaction mapping when ancillary carbon-bearing substrates may contribute
to the product. The 354 rows with missing precursor identities now also retain
537 carbon-bearing required-substrate alternatives from Terpedia; these are
explicit unresolved candidate structures, not assigned precursor identities.
RDKit mapped 282 of those alternatives completely and retained 14,252
candidate/inferred atom correspondences across the alternative layer; 255
alternatives and 3,629 product atoms remain unresolved.
An additional core-bridge projection joins 135 of these candidate precursor
structures to 25 exact-matched Terpedia core metabolite entities, spanning 28
CannabisDB target compounds. These gold dashed graph edges are a reviewable
bridge between the identity-set evidence and the core reaction inventory; they
remain outside balanced-reaction and CO₂-lineage metrics.
Each CannabisDB compound record carries its CO₂-lineage status and reachable
carbon-atom count, and NetworkDB links to the complete atom-level audit.
Every one of the 1,107 working reaction records also carries a carbon-mapping
summary from the RDKit report: 333 reactions are fully inferred, 281 retain
candidate mappings, 259 retain ambiguous carbon mappings, 234 retain
unresolved mappings, and 0 are unavailable. These statuses are
independent of enzyme status and are exposed as a separate Cytoscape filter.
Each reaction also records its inferred and candidate lineage-edge counts and
the source report used to derive them. The Pages Cytoscape view exposes
`non_enzymatic` as a separate evidence-status filter.
The generated `data/reports/carbon-mapping-work-queue.json` ranks the 493
reactions with unresolved or ambiguous product-carbon rows: 8,806 carbon
rows currently require mapping curation. Ranking prioritizes unresolved rows,
then ambiguous rows, then total blocked carbon weight. Each item retains
reaction source, mapping methods, known enzyme IDs, and attached genome
candidate proteins. Queue items now also retain reaction-direction evidence and
classify it as curated or raw-master orientation.
The lineage report also carries a direction-agnostic structural sensitivity
run. It is an upper bound for prioritizing reaction-direction curation, not a
physiological route: the directed result remains the authoritative pathway
view, and the reversible run never upgrades evidence status.
It also carries 2,264 candidate hypotheses and 4,063 deduplicated candidate
protein records; 946 reaction records have one or more attached candidate
proteins. All 4,063 candidate proteins are present in the 30,304-protein
UP000583929 reference-proteome FASTA and have sequence length/hash evidence in
`data/reports/genome-candidate-search.json`. This is proteome membership and
annotation evidence, not a homology score or confirmation of enzyme function.
The full proteome was also searched with DIAMOND blastp and specialized
sequence searches against reviewed reference proteins. 670 candidate
proteins currently have recorded hits and meet the screening threshold. These
thresholds define candidate evidence only and do not establish substrate
specificity or in-vivo activity.
The specialized IspD search for the added Rhea 13429 step found cannabis
protein A0A7J6EQA1 (gene G4B88_002492), with 63.4% identity, 77.9% query
coverage, and E-value 4.8e-97 to reviewed Arabidopsis IspD. This is a
sequence-supported candidate producer, not a confirmed cannabis enzyme.
The standalone `data/reports/enzyme-gap-audit.json` artifact audits the 12
enzyme-less reactions individually: all 12 have candidates, covering 27
candidate proteins; all 27 have reference-proteome sequence evidence, 7 have
DIAMOND hits, and 15 meet the current strong-candidate threshold. This audit
does not promote any candidate to a confirmed enzyme edge.

## Testable hypothesis set

`docs/data/testable-hypotheses.json` converts every candidate or blocked queue
record into a falsifiable hypothesis. Each record preserves the reaction and
exact participant IDs when available, candidate proteins and best sequence hit,
blocking causes, source provenance, and a proposed assay plan. The current
reaction set contains 2,264 records: 1,726 candidate records and 538 blocked
records; 1,347 are attached to a working reaction, including eight explicit
non-enzymatic conversion hypotheses. It also contains
6,220 target-level hypotheses, one for every CannabisDB compound, including
6,024 unresolved CO₂-lineage targets. TKS and OAC have dedicated coupled or
substrate-specificity assay plans in addition to recombinant-enzyme and
plant-validation steps. Target records retain CannabisDB identity fields and a
review priority: 88 explicit cannabinoid/Cannabis specialty records are high
priority, 152 exact-identity records are medium priority, and the remainder
stay normal priority until identity or route evidence improves. These are
future tests, not claims of confirmed activity. Each target hypothesis also
carries total, reachable, unresolved, and reversible-upper-bound carbon counts
plus the lineage blocker, so the assay queue can be ranked by missing carbon
rather than compound count alone.

## Terpedia forward hypothesis-edge layer

The GCP table `terpedia-489015.terpedia_core.terpene_hypothetical_forward_connections`
is published as `docs/data/hypothetical-forward-connections.json` and loaded as
a separate Cytoscape layer. It currently contributes 1,901 source-directed
connections (1,666 Rhea-derived and 235 MARTS-DB-derived) across 769 identity
compounds: 59 are bridged to CannabisDB by exact InChIKey and the remainder
are retained as Terpedia identity-set nodes. Purple dashed edges are candidate hypotheses
only; they are not included in balanced reaction, completeness, or CO₂ lineage
counts until source reaction, direction, atom balance, and enzyme evidence are
validated. Each record retains GCP source fields, structure-match mode, reaction
SMARTS, source enzyme/protein fields, and its claim boundary.
The complete companion reaction inventory contains 5,429 GCP hypotheses and
1,324 distinct products; products without a resolved corpus precursor remain
visible as explicitly unresolved product-inventory nodes.
The map also expands these records into 7,384 hypothesis-layer edges: 6,429
candidate edges and 955 red unresolved-substrate edges backed by 531 unique carbon-containing
missing-substrate structures. Exact canonical RDKit structure matching resolved
100 previously anonymous substrate structures, and unique connectivity matching
resolved a further 110 edge instances with explicit stereochemical uncertainty;
those edges remain candidate hypotheses and do not assert in-vivo activity.
Candidate-path atom reports also expose endpoint carbon deltas and the carbon
count of required or missing cosubstrates. This keeps methylation and other
multi-substrate steps unresolved when the new carbon is not itself connected
to the CO₂ lineage. The same report now includes a deduplicated
`carbon_source_gaps` queue: 22 candidate path steps are flagged, with 529
missing input-carbon atoms and one product-carbon deficit, for targeted
substrate/reaction curation.
The PubChem PUG View reconciliation report
(`docs/data/pubchem-chebi-xrefs.json`) covers 3,586 exact PubChem-resolved
CannabisDB records; 1,458 have explicit ChEBI cross-references, representing
1,854 ChEBI identifiers. It adds 15 new one-to-one candidate identity links
after existing exact, connectivity, XML-ChEBI, and name layers are respected.
These links remain candidate evidence and do not silently promote a structure
or pathway.
The GCP `terpene_enzyme_reaction_gene_evidence` table is preserved in
`docs/data/terpene-enzyme-reaction-gene-evidence.json` (4,639 records) and
attached to 514 hypothesis connections by MARTS reaction ID. It includes
source-species enzyme names, UniProt/GenBank identifiers, sequences, and assay
links; this evidence supports candidate annotation only and is not Cannabis
functional validation.
The companion `docs/data/terpene-biotransformation-enzyme-catalog.json`
preserves 585 curated Rhea enzyme-family records and joins 607 hypothesis-layer
edges by Rhea reaction ID, including EC, direction, family, and source links.
The companion `docs/data/terpedia-hypothesis-balance-audit.json` classifies
these SMARTS records as 1,762 balance-ready, 17 imbalanced, and 122 not
auditable; this is a review gate, not automatic promotion.
The companion `docs/data/terpedia-hypothesis-carbon-mapping.json` applies a
bounded RDKit exact-substructure correspondence to the substrate/product
identity structures: 388 connections are structurally inferred for 9,885
product carbon atoms, 12 retain candidate correspondences, and 1,501 remain
unresolved (39,621 product carbon atoms). Non-matches are retained explicitly;
the report does not use atom order or an unbounded MCS to fill gaps.

## Phase 1 balance gate

### Re-solved expanded candidate network (2026-09-04)

**Joint direction sensitivity:** disabling only the five reverse orientations
identified below, while keeping their source-written orientations and all
other candidate reactions available, yields **101** target certificates rather
than 108. None of the seven added targets retains a certificate; each is
solver-reported infeasible and no alternative was found. All 100 original
structure witnesses (101 target records) are preserved unchanged because none
uses a forbidden step. This is a model sensitivity result, not proof that the
reverse chemistry is impossible or that the seven metabolites are absent.

The audit reconsiders all 6,220 targets in the 1,588-equation candidate network.
There are 3,171 allowed directed steps after five are disabled, involving 1,667
distinct allowed reactant compounds. The target results are 101 certificates,
187 solver-reported infeasible, 5,897 without a producing candidate equation,
and 35 explicit exchange species. Target coverage and the reactant-compound
denominator are separate metrics. CO2 and the explicit carbon-free reservoir
are unchanged; no new organic seeds, thermodynamic constraints or atom traces
are introduced. Other reaction directions remain permissive hypotheses.

`data/reports/phase1-direction-sensitivity.json` records every target result,
the five precise forbidden step IDs and their review links, and references to
unchanged certificates in the expanded report. Terpedia stores 6,226 records in
`terpedia_core.cannabis_phase1_direction_sensitivity_20260904_v1`. This report
does not replace the unrestricted expanded map. Follow-up priorities are
direction-specific assays and alternative plant-supported reactions, not
further promotion of the unrestricted numerical count.

**Direction review:** all seven newly feasible target certificates use one or
more of the five newly selected candidate equations opposite to the source's
written equation. This is a source-orientation comparison through the explicit
published Rhea direction-family mapping, **not** a physiological direction
assignment or proof that the reverse reaction is impossible. The 108 numerical
certificates must not be read as 108 validated biosynthetic routes.

`data/reports/phase1-candidate-direction-review.json` preserves five complete
reaction reviews and seven target associations, with exact selected inputs,
outputs, extents, source equations, literature links, and candidate evidence
IDs. The new GCP table is
`terpedia_core.cannabis_phase1_candidate_direction_review_20260904_v1` (13
records including metadata). Each affected graph edge carries its review ID
and warning; the reaction panel lists discriminating tests.

The reviewed source masters are RHEA:25241 (ureidoglycine hydrolysis),
RHEA:27329 (NADH-linked urate oxidation), RHEA:16425 (formaldehyde oxidation),
RHEA:20001 (glycolaldehyde oxidation), and RHEA:64904 (GTP hydrolysis).
The last is used by the permissive model to synthesize GTP from GMP and
phosphate; energy feasibility is not established by element/charge balance.

The [HpxO primary paper](https://pubmed.ncbi.nlm.nih.gov/19260710/) characterizes
an FAD-dependent urate oxidase from *Klebsiella pneumoniae* in purine breakdown,
not a reverse Cannabis reaction. Its illustrated assays use NADPH; the exact
NADH/NADPH specificity of the proposed conversion needs separate checking.
The review retains the abstract/figure-caption evidence scope rather than
claiming a complete literature assessment.

Proposed tests compare each complete selected conversion with the source-written
direction, include inactive/no-protein controls, verify all products, and
discriminate cofactors and related substrates. Unsupported reverse activity
should prompt alternative plant-route curation, not an assertion that the
metabolite is absent. Original certificates and candidate evidence are retained
unchanged; only explicit risk annotations are added.

The [expanded candidate map](net.html?scenario=expanded) re-solves all 6,220
target records after adding the 116 catalog/backfill candidate equations to
the original 1,472-equation candidate network. Every admitted equation is
independently element-, isotope- and charge-balanced. No inferred completion
equations are added in this scenario. Candidate evidence is not proof of
Cannabis activity, specificity, physiological direction or complex assembly.

The 1,588-equation network has **108 target net-conversion certificates**
(107 exact structures), seven target records beyond the unchanged baseline
101. Newly feasible records are adenosine (CDB004791), guanosine (CDB004808),
inosine (CDB004818), uric acid (CDB004839), 5'-methylthioadenosine (CDB004887),
glycolaldehyde (CDB004953), and 5-amino-6-ribitylamino uracil (CDB004992).
Five newly candidate-linked equations occur in the selected certificates.

This differs from the **102** fully candidate-linked target certificates in
the frozen chemistry-only diagnostic: re-solving can select alternative
equations. These counts have different selection procedures and must not be
substituted for each other. Original baseline certificates remain verbatim.
The completion-sensitivity scenario is separate and unchanged.

CO2 remains the only carbon exchange; the same explicitly listed carbon-free
reservoir is allowed. Every certificate is checked with exact rational
arithmetic: no internal species is depleted overall, all positive net products
are retained, and net carbon input equals net carbon output. These are
stoichiometric hypotheses allowing regenerated pre-existing internal pools,
not proof of pool origin, energy feasibility, compartment compatibility, or
atom provenance. Both directions remain hypothetical.

Separately, all-reactants startup scope adds uric acid, for two nonseed target
records under the permissive carbon-free reservoir; 6,183 targets remain
startup-blocked and 35 are explicit exchange species. This qualitative scope
does not establish a quantitatively feasible startup sequence. The net audit
retains 5,897 targets with no producing candidate equation and 180 with
solver-reported infeasibility, not proof of biological absence.

The report is `data/reports/phase1-expanded-candidate-net.json`, with 7,021
export records for `terpedia_core.cannabis_phase1_expanded_candidate_net_20260904_v1`.
The static view is `docs/data/expanded-net-view/`; its manifest pins the report
and full evidence sources. Blue arrows identify newly screened candidate
equations; each projected arrow retains the complete reaction input/output
lists. The report and view preserve all target records, including gaps.

### Remaining catalog-gap source annotation audit (2026-09-04)

`data/reports/phase1-gap-annotations.json` audits all **349** equations still
missing candidate evidence after the combined catalog/backfill supplement.
It joins exact source reaction IDs through the pinned, published Rhea
direction-family table, then preserves complete source triples from 14
checksummed Rhea SPARQL responses for 341 master reactions.

- 200 equations have a source EC link; Cannabis function remains unresolved.
- 141 mapped equations have no source EC link; catalysis remains unresolved.
- 8 equations lack a published source-family mapping.
- 316 equations have source-linked literature citations available for review.

Missing EC annotations or reference proteins **do not establish spontaneous
chemistry**. This audit adds no enzyme, spontaneous-reaction, physiological
direction, or pathway-completion claims. Existing equations, full participant
lists, evidence and selected certificates remain unchanged. The ranking counts
membership in selected certificates, not reaction necessity.

The first-ranked gap, RHEA:46952 (EC 1.17.1.11), appears in 181 selected target
certificates. Its linked primary paper reports an electron-bifurcating formate
dehydrogenase complex purified from *Clostridium acidurici*, comprising four
subunits encoded by `hylCBA-fdhF2`; it does not demonstrate Cannabis activity.
See [Wang et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23872566/).
Follow-up must review organism, complex assembly, cofactors and measured
direction, and compare plant-supported alternative routes before treating this
as a high-priority missing Cannabis gene. A single homolog would not establish
the complete complex's catalytic capability.

The audit is exported as 705 provenance-preserving records for the versioned
Terpedia table `cannabis_phase1_gap_annotations_20260904_v1`.
Raw request metadata and responses are under `data/raw/phase1-gap-annotations/`;
the module `cannabis_carbon.phase1_gap_annotations` replays cached evidence.

### Plant purine annotation audit and protein hypotheses (2026-09-04)

`data/reports/phase1-plant-purine-references.json` preserves the complete
24-record response to the explicit Arabidopsis de novo purine UniProt query:
11 reviewed and 13 unreviewed annotations. This query is a pathway-discovery
aid, not proof of a complete pathway inventory. Its 13 annotated Rhea families
join exact balanced Terpedia equations: 10 already have candidate evidence,
and 3 lacked it. Q9SJ42 has no Rhea annotation and remains explicitly unjoined.
The source response, query, timestamps and hashes are retained for replay.

The three gaps were absent from earlier selected-certificate search queues.
This demonstrates a limitation of prioritizing only selected minimum-extent
routes: alternative plant-pathway reactions can remain unscreened even when
their equations are present in the catalog.

`data/reports/phase1-plant-purine-search.json` screens all **30,304 Cannabis
proteins** against all eight gap-reference sequences (all retrieved, all
unreviewed). Of 31 raw alignments, 16 pass the existing thresholds: at least
30% identity and 50% coverage of both query and reference. Four distinct
Cannabis proteins yield six protein–reaction hypotheses:

| Exact Rhea family | Reference-annotated activity | Cannabis candidates |
| --- | --- | --- |
| RHEA:16853 | Adenylosuccinate lyase | A0A7J6GQT2; A0A7J6I882 |
| RHEA:23920 | Adenylosuccinate lyase, second reaction | A0A7J6GQT2; A0A7J6I882 |
| RHEA:22192 | AICAR formyltransferase | A0A7J6E4V0; A0A7J6HKH3 |

Each hypothesis retains the exact balanced equation, source-family join,
reference review status, sequence provenance and passing alignment evidence.
Full participant lists are preserved. Multiple protein records are not counted
as distinct genes without locus evidence. Homology to an unreviewed reference
does not establish activity, substrate specificity, direction or compartment
in Cannabis. The next tests are exact-substrate product assays, controls for
the reference-annotated activities, and independent localization/cofactor
assessment; upstream precursor supply still needs pathway analysis.

These are separate discovery reports, not an automatic addition to an existing
candidate network. **No pathway-completeness count increases in this release.**
The expanded scenario remains at 108 conditional net target records, and its
five-direction sensitivity remains at 101; neither establishes physiological
flux or zero-pool startup. Atom tracing remains deferred.

The reports export 41 annotation-audit records and 33 screen records to
`terpedia_core.cannabis_phase1_plant_purine_references_20260904_v1` and
`terpedia_core.cannabis_phase1_plant_purine_search_20260904_v1`, respectively.
Raw public inputs are under `data/raw/phase1-plant-purine-references/` and
`data/raw/phase1-plant-purine-search/`. Replay modules have the corresponding
`cannabis_carbon.phase1_plant_purine_references` and
`cannabis_carbon.phase1_plant_purine_search` names.

### Exact purine-participant supply diagnostic (2026-09-04)

`data/reports/phase1-purine-precursor-audit.json` tests the **32 exact chemical
participants** of the 13 plant-annotated purine equations, plus the seven
previously added target structures. The 39 probes retain exact structures,
charges, coefficients and source reaction joins. All scenarios use the same
CO₂-only carbon exchange boundary and forbid the same five analyst-selected
reverse steps from the direction-sensitivity report.

| Scenario | Equations admitted | Participant net hypotheses | Explicit exchange participants | Solver-infeasible participants | Seven target net hypotheses |
| --- | ---: | ---: | ---: | ---: | ---: |
| Restricted expanded candidates | 1,588 | 6 | 5 | 21 | 0 |
| Same plus three unreviewed plant-reference hypotheses | 1,591 | 6 | 5 | 21 | 0 |
| Restricted full balanced catalog, regardless of enzyme evidence | 13,995 | 27 | 5 | 0 | 7 |

This local comparison points to missing candidate-linked chemistry in the
current model: the three newly screened purine reactions alone do not restore
precursor supply, whereas additional catalog chemistry permits net conversion.
It does **not** establish that all catalog routes operate in Cannabis, that
missing genes exist, or that these reactions are the only possible alternatives.
Other reaction directions remain hypothetical and internal pools may be
pre-existing and regenerated. Each probe is solved separately; this is not a
joint-supply, pathway-ordering, energetic-feasibility or startup certificate.

The report retains every successful exact net certificate, all 152 selected
or plant-annotated equations, and 230 participating or exchange structures.
Its **70-reaction candidate-evidence queue** records full equations and exact
selected uses in catalog certificates for probes unresolved in the augmented
candidate model. Queue membership means no candidate link in this model, not
absence of a protein or unique reaction necessity. Before another genome
search, audit existing search reports and source annotations for these exact
equations, including direction, reference coverage and complex requirements.

The 459-record export is stored in Terpedia as
`terpedia_core.cannabis_phase1_purine_precursor_audit_20260904_v1`.
`cannabis_carbon.phase1_purine_precursor_audit` reproduces the diagnostic;
`add_gap_queue` derives the queue from saved certificates without re-solving.
Tests independently replay every selected equation and successful certificate,
including element/isotope/charge balance, exact net coefficients, permitted
exchanges and forbidden directions. Published whole-metabolome completeness
counts and map scenarios remain unchanged. Atom tracing remains deferred.

### Previously unscreened purine-alternative reactions (2026-09-04)

`data/reports/phase1-purine-gap-references.json` audits all 70 reactions in the
restricted catalog precursor queue against nine pinned prior search reports.
It preserves previous results for 48 equations and identifies **22 previously
unscreened equations**, all with explicit published Rhea family mappings.
One successful reviewed, nonfragment UniProt query returns 833 exact-family
reference proteins covering 16 equations. Six equations return no reviewed
reference; that is an annotation gap, not evidence of spontaneous chemistry or
biological absence. Raw annotations, query provenance and checksums are retained.

`data/reports/phase1-purine-gap-search.json` screens all **30,304 Cannabis
proteins** against all 833 retrieved reference sequences. Fourteen sequence
retrieval batches completed without failure. Of 1,247 raw alignments, 257 pass
the existing 30% identity and 50% query/reference coverage thresholds.

| Equation search outcome | Equations |
| --- | ---: |
| Candidate proteins found | 10 |
| Weak hits only | 4 |
| No hits despite available references | 2 |
| No reference sequence | 6 |

The passing alignments identify **29 distinct Cannabis proteins** and **36
protein–reaction hypotheses** across the ten equations. Four focused targets
have a selected catalog certificate containing at least one of these new
candidate-linked reactions; this is not a count of rescued pathways.
Each equation retains its full participant lists and exact source-family joins;
each candidate retains sequence and alignment provenance. Review substrate
specificity, catalytic residues, domains, complex partners and direction before
designing exact-input/product assays with reference and no-enzyme controls.

This release adds evidence reports, not an automatic change to a pathway
scenario. Net conversion must be re-solved after admitting the new hypotheses;
the full-catalog certificates are not automatically candidate-supported.
The pathway-inference evidence rules keep reviewed reference annotation,
Cannabis homology, biochemical activity and pathway execution separate.
Published map completeness counts remain unchanged. Atom tracing is deferred.

Terpedia tables `cannabis_phase1_purine_gap_references_20260904_v1` (927 records)
and `cannabis_phase1_purine_gap_search_20260904_v1` (1,156 records) preserve these
artifacts in `terpedia_core`. Replay modules are
`cannabis_carbon.phase1_purine_gap_references` and
`cannabis_carbon.phase1_purine_gap_search`; raw public inputs use the corresponding
`data/raw/phase1-purine-gap-references/` and `data/raw/phase1-purine-gap-search/`
directories. Prior results are retained rather than overwritten.

### Combined purine candidates: whole-metabolome recalculation (2026-09-05)

`data/reports/phase1-purine-candidate-net.json` adds all 13 newly screened
equations to the prior 1,588-equation candidate network. The 1,601 admitted
equations are independently element/isotope/charge balanced. The supplement
contains 33 distinct Cannabis proteins, with three unreviewed plant-reference
equations kept distinct from ten reviewed-reference equations. Neither class
establishes Cannabis enzyme activity.

Both scenarios retain and assess **all 6,220 target records** under the same
explicit CO₂-only carbon exchange boundary. Existing exact certificates are
preserved verbatim only when they obey the scenario's direction restrictions;
all other target structures are re-solved.

| Scenario | Conditional net target records | Exact structures | Solver-infeasible records | No net-producing candidate equation | Explicit exchange records |
| --- | ---: | ---: | ---: | ---: | ---: |
| All candidate directions hypothetical | 109 | 108 | 179 | 5,897 | 35 |
| Five reviewed reverse steps forbidden | 101 | 100 | 187 | 5,897 | 35 |

The single new permissive target is **5′-deoxyadenosine, CDB004932**. Its selected
certificate uses four newly admitted equations, including two unreviewed
plant-reference hypotheses. It also uses the previously flagged reverse steps
for RHEA:25241 and RHEA:27329. The five-step restriction removes this new route
and the previous seven permissive additions; no replacement target certificates
are found in the restricted candidate model. Source equation ordering alone is
not a physiological direction constraint. These results do not establish
absence of alternative biology, enzyme function, feasible energy coupling or
zero-pool startup. Internal pools may be pre-existing and regenerated.

The static map now exposes the [purine supplement](net.html?scenario=purine)
and [five-step restricted supplement](net.html?scenario=purine-restricted), using
one shared data bundle. The permissive view opens on 5′-deoxyadenosine by default.
Full equations, source links and candidate evidence remain attached to arrows;
the new unreviewed-reference highlight and existing direction-risk highlight
dim other arrows without deleting participants. Startup is explicitly marked
as not recomputed in these views. Earlier baseline and expanded scenarios remain
available unchanged.

The report's 13,295-record export is stored as
`terpedia_core.cannabis_phase1_purine_candidate_net_20260905_v1`.
Replay modules are `cannabis_carbon.phase1_purine_candidate_net` and
`cannabis_carbon.phase1_purine_candidate_view`. Tests independently check all
selected equations and certificates, original-target preservation, source hashes,
reference confidence, both graph projections and scenario switching. Net
certificates are conditional hypotheses, not confirmed pathway completeness or
atom provenance. Atom tracing remains deferred.

### Recovering skipped reference discovery (2026-09-05)

`data/reports/phase1-deferred-references.json` audits all **60** selected
purine-alternative gaps still absent from the current 1,601-equation candidate
network. It joins each prior search row to its checksummed discovery record,
rather than equating presence in a search report with a completed search.
For **seven equations**, every prior attempt had `no-reference-sequence` because
discovery was explicitly `not-searched-in-priority-pass`. Those are workflow
omissions, not negative reference searches. The other 53 gaps retain their
existing evidence; failed lookups, genuine empty responses and weak/no-hit
results are not reclassified as skipped.

The missing reviewed-reference query found 29 proteins for five equations;
two equations returned no reviewed references. The complete-proteome screen in
`data/reports/phase1-deferred-search.json` retrieved all 29 references and
searched all 30,304 Cannabis proteins. It records 266 raw alignments and 67
passing alignments, with the following equation outcomes:

- 1 equation with nine candidate proteins;
- 2 equations with weak hits only;
- 2 equations with no hits despite available references;
- 2 equations with no reference sequence.

The candidate equation has source variants RHEA:10181 and RHEA:10182, joined
through the published family for 10-formyltetrahydrofolate dehydrogenase.
**Full catalytic architecture remains unestablished.** The nine Cannabis
sequences are 466–543 residues long; each highest-bitscore representative match
is to a 902-residue reference and covers only 52.4–54.0% of that reference.
These are partial-reference alignment leads, not evidence for an intact
multidomain catalytic system. Alignment coordinates/domain localization were
not established by this screen. Review complete domain composition, required
partners and exact reaction direction before admitting the reaction to a
candidate pathway or designing an activity assay.

The three highest-ranked formerly skipped gaps remain unresolved: RHEA:33871
returned no reviewed reference, while RHEA:19029 and RHEA:20896 produced only
weak hits. Each occurs in 20 selected precursor certificates, which is a
prioritization measure, not evidence of reaction necessity.

The 98-record reference/audit export and 114-record screen export are stored as
`terpedia_core.cannabis_phase1_deferred_references_20260905_v1` and
`terpedia_core.cannabis_phase1_deferred_search_20260905_v1`. Replay modules are
`cannabis_carbon.phase1_deferred_references` and
`cannabis_carbon.phase1_deferred_search`. Historical reports remain unchanged;
tests require explicit skipped-discovery evidence before recovering a row.
This release changes neither map scenarios nor pathway counts. Atom tracing
remains deferred.

### Ureidoglycine direction review and alternative-route sensitivity

The reviewed RHEA:33871 equation describes ureidoglycine decay, but all 20
selected precursor certificates that used it required the opposite direction.
The [primary paper abstract](https://www.nature.com/articles/nchembio.445)
reports spontaneous decay supplying glyoxylate in a bacterial transamination
study; it does not establish reverse synthesis or Cannabis activity. Missing
EC or reviewed sequence alone is not evidence of spontaneous chemistry.
The exact identities, source-direction joins, literature-access limitations
and proposed experiments are preserved in
`data/curation/ureidoglycine-direction-review.json`.

`phase1-decay-sensitivity.json` independently re-solves all 39 exact probes
from the precursor audit in the 13,995-equation balanced catalog. It preserves
the five previous direction restrictions and excludes only reverse
ureidoglycine synthesis as a sixth sensitivity constraint; decay remains
available. CO2 remains the sole carbon exchange.

All 20 affected probes have alternative exact net-conversion certificates.
Overall, 34 probes have certificates and five are explicit exchange species;
all seven focused target structures retain chemistry-only certificates.
The selected routes use 139 equations and contain 66 equations without
candidate links in the current 1,601-equation model. These are witness-specific
gaps, not required plant reactions or a change in whole-metabolome completeness.
Other directions, regenerated internal pools, compartments and enzyme activity
remain unverified. No candidate is promoted by this sensitivity analysis.

The 474-record export is stored in
`terpedia_core.cannabis_phase1_decay_sensitivity_20260905_v1`.
Replay: `python -m cannabis_carbon.phase1_decay_sensitivity`.
Tests independently check exact stoichiometry, isotope/element/charge balance,
CO2-only carbon input, all six restrictions, source hashes and gap membership.
Existing map scenarios and atom-accounting artifacts remain unchanged.

### Protein searches for newly exposed replacement-route gaps

The replacement-route inventory contains five equations without prior search
records and 61 equations with prior evidence. All 61 retain their original
results, including the withheld partial-reference RHEA:10180 leads; no old
negative or withheld candidate is silently promoted.

Reviewed, nonfragment UniProt queries returned 524 references for all five new
equations. All sequences were retrieved and screened against all 30,304
Cannabis proteins. The screen produced 1,068 alignments, of which 434 pass the
existing identity/coverage/e-value thresholds. Four Cannabis proteins provide
four reaction hypotheses covering two equations:

- RHEA:16845: A0A7J6ER78 and A0A7J6ETR9 match a reviewed plant citrate-synthase
  reference (P49298), with 87.9–88.2% identity and 100% reference coverage.
  **The selected replacement witness uses the equation opposite to its
  source-written synthesis direction.** Homology does not validate this
  reverse step or establish a citrate-cleavage mechanism.
- RHEA:20780: A0A7J6EDY8 and A0A7J6HES8 match reviewed Arabidopsis guanylate
  kinases (Q9M682 and Q94JM2). Their highest-bit-score alignments have
  65.0–66.6% identity and 93.3–97.7% reference coverage. Localization and
  Cannabis activity remain untested; the source-written direction matches
  this selected witness but is not proof of physiological flux.

RHEA:35799, RHEA:25860 and RHEA:24836 returned no hits under this screen.
This is not evidence of biological absence. Neither new candidate equation
has been admitted to the 1,601-equation candidate network, and no pathway or
completeness count changes in this release. These are sequence-supported
test leads, not characterized enzymes or distinct-gene assertions.

Reports and replay modules are `phase1-replacement-references` and
`phase1-replacement-search` (`cannabis_carbon` module namespace).
The 597-record discovery/audit export and 977-record search export are stored
in `terpedia_core.cannabis_phase1_replacement_references_20260905_v1` and
`terpedia_core.cannabis_phase1_replacement_search_20260905_v1`.
Raw requests, reference sequences, alignments, sequence checksums and exact
source-reaction joins are retained for replay. Atom tracing remains deferred.

### Whole-metabolome test of the replacement candidates

`phase1-replacement-candidate-net.json` evaluates the two new sequence-supported
equations as additions to the 1,601-equation candidate model. The resulting
1,603-equation test model retains all prior evidence IDs and adds four protein
hypotheses without upgrading their evidence class to characterized activity.
All 6,220 target records are accounted for in both scenarios:

- Permissive directions: 109 target records (108 exact structures) retain
  net-conversion certificates; 179 are solver-infeasible, 5,897 have no
  net-producing candidate equation, and 35 are explicit exchange species.
- Six reverse steps forbidden: 101 target records (100 exact structures)
  retain certificates; 187 are solver-infeasible, with the same 5,897
  no-producer records and 35 exchange species. This preserves the five prior
  restrictions and additionally excludes reverse citrate-synthase use.

Neither scenario gains a target certificate. Existing valid certificates are
preserved exactly; every other exact target structure is reconsidered with the
augmented model. No internal depletion or organic carbon import is accepted.
This is an enzyme-candidate inventory gain, not improved pathway completeness.
Reverse ureidoglycine synthesis is not an extra constraint here because its
equation remains absent from this candidate model altogether.

The report contains both new reaction evidence records even though neither
occurs in a preserved certificate. It does not replace the historical static
map scenarios or relabel their 1,601-equation inventory as 1,603. The 13,267-row
export is stored in
`terpedia_core.cannabis_phase1_replacement_candidate_net_20260905_v1`.
Replay: `python -m cannabis_carbon.phase1_replacement_candidate_net`.
Tests independently balance all 1,603 candidate equations and replay every
saved certificate using exact rational stoichiometry. Startup, physiological
direction, compartmentation and flux remain unestablished; atom tracing is
still deferred.

### Whole-inventory no-producer and identity audit

`phase1-no-producer-audit.json` retains all 5,897 no-candidate-producer target
records (5,883 exact encoded structures) from the 1,603-equation model.
The categories below use those current structures and do not resolve source
disagreements:

| Gap category | Target records | Next action |
| --- | ---: | --- |
| Exact full-catalog producer, no admitted candidate enzyme | 120 | Review 208 distinct producing equations and enzyme evidence |
| No exact participation, but diagnostic structure-variant leads | 588 | Resolve identity/protonation/stereochemistry without merging |
| No exact or diagnostic catalog match | 5,188 | Find or infer a complete balanced producing reaction |
| Exact catalog participation only with zero net change | 1 | Do not mistake participation for production |

The zero-net-only record is CDB000139, Aldehydo-L-rhamnose. Full-catalog
production is an upper bound considering both hypothetical directions, not
proof of upstream supply, Cannabis activity or a CO2 route.

Diagnostic keys separately remove stereochemistry, apply RDKit Uncharger, or
combine the two, retaining isotope labels. They yield leads for 139, 395 and
609 target records respectively; these overlapping counts must not be added.
All matching structures are retained, without picking a preferred identity or
creating reaction edges. Uncharger-key equality is not a validated proton
transfer, pH model, tautomer equivalence or enzyme annotation.

The audit separately compares all 6,220 SDF-derived and XML accessions and
retains both assertions for all 71 structure disagreements. Identity resolution
takes precedence over pathway inference for the 64 affected no-producer
records. XML structure/name/formula/InChIKey fields are checked against the
archived XML rather than reconstructed from names. Original xrefs remain
source claims and may be internally inconsistent.

The 17,356-record export is stored in
`terpedia_core.cannabis_phase1_no_producer_audit_20260905_v1`.
Replay: `python -m cannabis_carbon.phase1_no_producer_audit`.
No reaction, structure, source assertion or completeness metric is promoted
by this audit. Atom tracing remains deferred.

### Registry review of the 71 source identity conflicts

`phase1-identity-conflict-review.json` preserves 142 separate source assertions
and queries both structures by freshly computed InChIKey. Live Terpedia
`terpene_identity_set` corroborates 14 assertions and PubChem corroborates 132;
these are assertion counts, not independent compounds or validated Cannabis
occurrence. All nine lookups completed. Nine source-reported InChIKeys disagree
with their own SMILES and must not silently drive identity joins.

Of the 71 pairs, 48 share a standardized InChIKey despite distinct encoded
structures, 10 share only the connectivity key, and 13 differ in connectivity.
Standardized-key equality is explicitly separated from exact encoded-structure
agreement. It does not authorize merging tautomer, stereochemical or other
representation distinctions.

Four priority name queries yield these provisional findings:

| Record | Registry-supported finding | Migration status |
| --- | --- | --- |
| CDB006156 Glycerol | XML structure exactly matches [PubChem 753](https://pubchem.ncbi.nlm.nih.gov/compound/753); SDF structure exactly matches [(+)-beta-irone, 10219919](https://pubchem.ncbi.nlm.nih.gov/compound/10219919) | XML structure provisionally supported for this name; no historical overwrite |
| CDB000142 D-arabitol | SDF structure exactly matches the [D-arabitol name result, 94154](https://pubchem.ncbi.nlm.nih.gov/compound/94154) | SDF structure provisionally supported; XML alternative retained |
| CDB000546 Acetamide | SDF structure exactly matches [PubChem 178](https://pubchem.ncbi.nlm.nih.gov/compound/178); XML form shares its standard key but not exact encoded structure | Keep distinct forms; no implicit tautomer reaction |
| CDB006169 Ribitol | Neither source structure exactly matches the [name-query structure, 6912](https://pubchem.ncbi.nlm.nih.gov/compound/6912); XML is a different-connectivity hydrocarbon | Unresolved; stereochemical review required |

These are registry-supported proposals for versioned reconciliation, not
primary Cannabis-source corrections, biological observations or new pathways.
The three existing conditional certificates remain attached to their original
stored structures. Both copies of the source external IDs originate from XML
enrichment and are not independent SDF/XML corroboration. No graph structure,
carbon accounting or completeness metric is changed in this release.

The 227-record export is stored in
`terpedia_core.cannabis_phase1_identity_conflict_review_20260905_v1`.
Replay module: `cannabis_carbon.phase1_identity_conflict_review`.
Raw Terpedia SQL and PubChem requests/responses are checksum-pinned, including
negative exact-key results. Field consistency, registry matches, exact
structure comparisons and unresolved alternatives remain separately inspectable.

### Versioned identity branches: priority source conflicts

`data/reports/phase1-identity-branches.json` retains both source assertions for
four priority accessions and tests all eight exact structures independently in
both the 1,603-equation permissive and six-reverse-step-restricted models.
No historical target, structure, carbon inventory or completeness count is replaced.
Three field-level named-structure choices remain provisional; neither source's
names, formulas or external IDs are inherited wholesale.

| Accession | Source structure branch | Result in both models |
| --- | --- | --- |
| CDB006156, Glycerol | XML, C3H8O3; provisional named structure, PubChem 753 | Exact net CO₂-conversion hypothesis; 56 catalog producing equations |
| CDB006156, Glycerol | SDF-derived C14H22O structure | No candidate producing equation; no exact catalog producer |
| CDB000142, D-arabitol | SDF-derived provisional named structure; XML kept separately | Both exact structures have net hypotheses; five catalog producing equations each |
| CDB000546, Acetamide | SDF-derived provisional named structure | Exact net hypothesis; two catalog producing equations |
| CDB000546, Acetamide | XML encoded tautomer | No exact structure in the balanced network; no implicit tautomer conversion added |
| CDB006169, Ribitol | SDF-derived structure, identity unresolved | Exact net hypothesis; five catalog producing equations, **not proof of ribitol production** |
| CDB006169, Ribitol | XML hydrocarbon structure, identity unresolved | No candidate producing equation; no exact catalog producer |

The 16 branch/scenario tests produce ten exact net certificates, four
no-candidate-producer results and two no-exact-network-structure results.
These are alternative-structure tests, **not ten newly covered metabolites**.
All certificates retain exact rational reaction extents, imports, exports,
regenerated pools and reaction-level enzyme-evidence IDs. Only CO₂ supplies net
carbon. Regenerated pre-existing pools, hypothetical directions, carbon-free
exchanges, physiological flux and Cannabis occurrence remain explicit limitations.
Catalog producer counts use either chemical orientation, not established biological
direction or all-input availability. Atom tracing remains deferred.

Replay with `PYTHONPATH=src ./.venv/bin/python -m cannabis_carbon.phase1_identity_branches`.
The versioned Terpedia export is `data/derived/phase1-identity-branches.ndjson`;
all 65 records were read back and matched exactly in
`terpedia_core.cannabis_phase1_identity_branches_20260905_v1`
(`data/reports/phase1-identity-branches-gcp.json`).
the report records checksums of its registry review, balanced network and candidate
model inputs. Before promoting a provisional structure, verify the primary Cannabis
identification evidence; ribitol additionally requires stereochemical reconciliation.

### Working-network balance audit

The Phase 1 audit is stored in
`data/reports/phase1-balance-audit.json`:

- 727 of 1,107 working reactions are fully element- and charge-balanced after
  the RDKit canonical-structure and concrete reaction-SMILES fallbacks are applied.
- 0 reactions are explicitly imbalanced after computation.
- 380 reactions remain not auditable because their participants lack sufficient
  formula/structure or charge information.
- The source-only fields report 676 element-balanced and 801 charge-balanced
  reactions; computed fields are retained separately for provenance.

Phase 1 accepts a reaction only as stoichiometrically balanced when both element
and charge checks pass. When source formulas are absent, the audit may derive
formula and formal charge from Terpedia canonical SMILES or a concrete Rhea
reaction SMILES with RDKit; fallback results remain distinct from source
assertions. This gate does not claim enzyme function, reaction direction,
carbon atom provenance, or biological flux.
Rhea generic compounds containing wildcard substituents are classified as not
auditable rather than being assigned guessed element counts.

The companion integrity report `data/reports/artifact-validation.json` is
generated by `cannabis-carbon validate-artifacts`. It currently passes for all
6,220 CannabisDB compounds (285,623 carbon atoms), all 1,107 reaction mapping
rows, all atom evidence-field checks, and the no-imbalanced-reactions check.
This is an artifact-completeness gate, not evidence that every reaction is
biologically active.
