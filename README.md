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

CannabisDB currently provides bulk downloads rather than a public API:

`https://cannabisdatabase.ca/simple/download_compound_as_sdf`

`https://cannabisdatabase.ca/simple/download_compound_as_xml`

`https://cannabisdatabase.ca/simple/download_protein_as_xml`

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
