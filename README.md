# Cannabis carbon provenance

This repository is an attempt to make the cannabis metabolome carbon-complete:
for every known cannabis metabolite, identify a defensible path from inorganic
carbon (`CO2`) to each carbon atom.

The project is deliberately evidence-first. A structure match is not a pathway,
and an enzyme annotation is not proof of flux. Every inferred edge must retain
its source, reaction direction, atom-mapping method, and uncertainty.

## Current scope

- Cannabis Compound Database (CannabisDB) is the metabolite/protein source.
- RDKit performs molecular parsing, formula audits, and atom-level mapping.
- Reactions are represented as mapped reactant/product structures with explicit
  carbon provenance. Unmapped product carbons are reported as gaps.
- Enzyme evidence is kept separate from structural atom mapping.
- Missing reaction enzymes trigger genome-wide candidate discovery against
  annotated cannabis proteins, using enzyme-family homology, catalytic motifs,
  domain architecture, localization, and expression where available.
- `docs/data/networkdb.json` is the unified NetworkDB snapshot: all 6,220
  CannabisDB compounds, all 1,253 working Terpedia metabolite nodes, all 1,093
  working reactions, participant coefficients, enzyme associations, exact/candidate
  identity links, and reaction-level carbon-mapping summaries are retained in
  one source-linked artifact.
- `docs/data/testable-hypotheses.json` turns the unresolved queue into 2,243
  falsifiable reaction records and adds one target-level hypothesis for each
  of the 6,220 CannabisDB compounds, with candidate proteins, required inputs,
  blockers, provenance, and proposed recombinant/plant validation assays. The
  target records retain CannabisDB names, aliases, formulas, SMILES, and review
  priority so future experiments can be selected without losing identity context.

## Genome-to-enzyme discovery

For each carbon-producing reaction without a known cannabis enzyme, the planned
evidence ladder is:

1. Resolve the reaction and exact substrate/product structures from Terpedia.
2. Find characterized enzymes for the same EC/Rhea chemistry in curated public
   sources and build a reference sequence set.
3. Search cannabis protein translations from a pinned genome annotation (the
   Ensembl Plants `cs10` assembly and NCBI RefSeq are supported sources).
4. Filter and rank candidates by profile homology, catalytic residues/motifs,
   domain completeness, subcellular targeting, trichome/tissue expression, and
   genomic neighborhood.
5. Test each candidate against the RDKit atom-mapping requirements of the
   reaction. Sequence similarity alone cannot establish substrate specificity
   or in-vivo flux.

Candidate results will be labeled `homology_candidate`, `annotation_supported`,
`biochemically_supported`, or `rejected`, with the first two never promoted to
confirmed pathway edges automatically.

## Runtime split

The Cytoscape visualization is a static GitHub Pages artifact and does not
require GCR, Docker, or Cloud Run. GCP is reserved for data ingestion, RDKit
mapping, genome searches, and report generation; those jobs publish versioned
JSON artifacts consumed by Pages.

CannabisDB currently provides bulk downloads rather than a public API:

`https://cannabisdatabase.ca/simple/download_compound_as_sdf`

`https://cannabisdatabase.ca/simple/download_compound_as_xml`

`https://cannabisdatabase.ca/simple/download_protein_as_xml`

Terpedia sources and known gaps are tracked in
[`data/terpedia-sources.json`](data/terpedia-sources.json) and
[`data/reaction-gaps.json`](data/reaction-gaps.json). The current Terpedia KB
contains reaction-oriented tooling, but the public KB API was not reachable in
the deployment environment; the next ingestion step must export or expose its
RHEA/EC reaction records to this project.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
python -m cannabis_carbon download --out data/raw
python -m cannabis_carbon inspect-sdf data/raw/compounds.sdf
```

The download command records the retrieval timestamp and SHA-256 checksums.
Raw database files are intentionally gitignored; generated reports are
versionable and should include the source manifest.

## Status vocabulary

`direct`: atom mapping is supported by an explicitly mapped reaction.

`inferred`: mapping is proposed by RDKit under a stated transformation rule.

`unresolved`: one or more product carbons have no defensible precursor carbon.

`not-a-pathway`: a database association (for example, a protein target) does
not establish biosynthetic production.

## Scientific boundary

No finite experiment can establish a universal cannabis metabolome across all
cultivars, tissues, developmental stages, environments, and processing states.
This project therefore aims for a versioned, auditable carbon-provenance graph,
with confirmed, candidate, and unresolved branches kept distinct.
