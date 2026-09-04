# Terpedia data reference

This file is the source contract for the cannabis carbon-provenance project.
It records what has actually been imported, how entities are identified, and
which fields are required before a reaction can contribute carbon provenance.

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
