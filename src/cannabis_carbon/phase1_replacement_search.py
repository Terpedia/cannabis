"""Screen new replacement-route references against the whole Cannabis proteome."""
from pathlib import Path
from .phase1_new_protein_search import run as search, export_table


def run():
    output = Path('data/reports/phase1-replacement-search.json')
    search(source=Path('data/reports/phase1-replacement-references.json'),
           raw=Path('data/raw/phase1-replacement-search'), output=output,
           claim_boundary='Whole Cannabis proteome search for previously unsearched equations in replacement chemistry witnesses. Homology does not establish exact activity, direction, complex assembly, compartment or CO2 pathway execution. Prior partial-reference leads remain withheld and historical scenarios unchanged. Atom tracing remains deferred.')
    export_table(output, Path('data/derived/phase1-replacement-search.ndjson'))


if __name__ == '__main__':
    run()
