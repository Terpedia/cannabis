"""Separate unreviewed-plant-reference candidate screen; never promote activity."""
from pathlib import Path
from .phase1_new_protein_search import run as search, export_table


def run():
    source=Path('data/reports/phase1-plant-purine-references.json')
    output=Path('data/reports/phase1-plant-purine-search.json')
    search(source=source,raw=Path('data/raw/phase1-plant-purine-search'),output=output,
        evidence_class='unreviewed-Arabidopsis-purine-reference-homology-candidate',
        additional_blockers=('unreviewed-reference-annotation','plant-reference-does-not-establish-Cannabis-function'),
        claim_boundary='Whole Cannabis proteome screen of explicitly joined plant purine-pathway reference annotations. Gap references are unreviewed and retained as such. Homology does not establish exact activity, direction, compartment or a CO2 pathway. Existing candidate networks and certificates remain unchanged. Atom tracing is deferred.')
    export_table(output,Path('data/derived/phase1-plant-purine-search.ndjson'))


if __name__=='__main__': run()
