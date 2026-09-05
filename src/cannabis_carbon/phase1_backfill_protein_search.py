"""Screen references recovered from previously skipped reaction-family lookups."""
from pathlib import Path
from .phase1_new_protein_search import run as search, export_table


def run():
    result = search(source=Path('data/reports/phase1-reference-backfill.json'),
        raw=Path('data/raw/phase1-backfill-protein-search'),
        output=Path('data/reports/phase1-backfill-protein-search.json'),
        evidence_class='reference-backfill-direction-unresolved-homology-candidate',
        additional_blockers=('chemistry-only-net-certificate-not-physiological-pathway',),
        claim_boundary='Full Cannabis proteome screen of reviewed reference leads recovered by completing previously skipped reaction-family lookups. Prior no-reference-sequence status was not proof of a negative reference search. Homology is not experimental activity, exact specificity, physiological direction or CO2 pathway confirmation. Existing screens, graph evidence and net chemistry are unchanged. Atom tracing remains deferred.')
    export_table(Path('data/reports/phase1-backfill-protein-search.json'),
        Path('data/derived/phase1-backfill-protein-search.ndjson'))
    return result


if __name__ == '__main__':
    run()
