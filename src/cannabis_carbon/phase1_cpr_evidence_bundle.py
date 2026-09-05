"""Lossless reductase candidates and explicit FNS-II carrier-interface gaps."""
from .phase1_reference_gap_bundle import run

SOURCES = (
    'data/reports/phase1-cpr-annotation-audit.json',
    'data/curation/cpr-fnsii-carrier-interface-review.json',
    'data/reports/phase1-cpr-references.json',
    'data/reports/phase1-cpr-search.json',
    'data/reports/phase1-cpr-annotations.json',
)


if __name__ == '__main__':
    run(SOURCES, 'cpr-evidence-bundle')
