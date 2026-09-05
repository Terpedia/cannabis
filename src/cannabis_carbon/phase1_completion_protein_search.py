"""Whole-Cannabis-proteome screen using original MARTS, not template, references."""
from pathlib import Path
from .phase1_completion_protein_discovery import BOUNDARY
from .phase1_new_protein_search import run


if __name__ == '__main__':
    run(source=Path('data/reports/phase1-completion-protein-discovery.json'),
        raw=Path('data/raw/phase1-completion-protein-search'),
        output=Path('data/reports/phase1-completion-protein-search.json'),
        evidence_class='original-MARTS-source-homology-for-inferred-stoichiometry',
        additional_blockers=('inferred-inorganic-stoichiometry-unverified',
                             'original-MARTS-exact-product-identity-unverified'),
        claim_boundary=BOUNDARY)
