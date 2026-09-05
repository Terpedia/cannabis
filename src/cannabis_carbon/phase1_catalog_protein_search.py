"""Screen newly discovered catalog-gap references against all Cannabis proteins."""
from pathlib import Path

from .phase1_new_protein_search import run as search, export_table


def run():
    result = search(
        source=Path('data/reports/phase1-catalog-references.json'),
        raw=Path('data/raw/phase1-catalog-protein-search'),
        output=Path('data/reports/phase1-catalog-protein-search.json'),
        evidence_class='catalog-net-gap-direction-unresolved-reference-homology',
        additional_blockers=('chemistry-only-certificate-not-physiological-pathway',
            'internal-pool-origin-unestablished'),
        claim_boundary='Screen of exact published reaction-family reference leads for previously unscreened chemistry-only catalog net gaps against the full Cannabis reference proteome. Homology is not activity, exact substrate specificity, physiological direction, or an established CO2 pathway. Prior screens and baseline certificates are unchanged. Atom tracing remains deferred.')
    export_table(Path('data/reports/phase1-catalog-protein-search.json'),
        Path('data/derived/phase1-catalog-protein-search.ndjson'))
    return result


if __name__ == '__main__':
    run()
