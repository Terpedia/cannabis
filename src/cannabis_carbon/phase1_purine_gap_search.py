"""Whole-proteome screen for previously unscreened purine-alternative gaps."""
from pathlib import Path
from .phase1_new_protein_search import run as search, export_table


def run():
    output = Path('data/reports/phase1-purine-gap-search.json')
    search(source=Path('data/reports/phase1-purine-gap-references.json'),
           raw=Path('data/raw/phase1-purine-gap-search'), output=output,
           claim_boundary='Reviewed exact reaction-family reference homology screen for previously unscreened gaps in selected restricted chemistry-only purine certificates. All 30,304 Cannabis proteins are searched. No hit or absent reference does not establish biological absence. Candidates do not establish activity, exact specificity, direction, joint precursor supply or CO2 pathway completion. Existing map scenarios and certificates are unchanged. Atom tracing remains deferred.')
    export_table(output, Path('data/derived/phase1-purine-gap-search.ndjson'))


if __name__ == '__main__':
    run()
