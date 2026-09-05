"""Screen references recovered from explicitly skipped discovery steps."""
from pathlib import Path
from .phase1_new_protein_search import run as search, export_table


def run():
    output = Path('data/reports/phase1-deferred-search.json')
    search(source=Path('data/reports/phase1-deferred-references.json'),
           raw=Path('data/raw/phase1-deferred-search'), output=output,
           claim_boundary='Whole Cannabis proteome screen for seven equations whose previous reference discovery was explicitly skipped. Historical records are preserved, not recast as biological negatives. Reviewed reference homology does not establish Cannabis enzyme activity, specificity, direction, complex assembly or CO2 route execution. Published candidate network scenarios remain unchanged. Atom tracing remains deferred.')
    export_table(output, Path('data/derived/phase1-deferred-search.ndjson'))


if __name__ == '__main__':
    run()
