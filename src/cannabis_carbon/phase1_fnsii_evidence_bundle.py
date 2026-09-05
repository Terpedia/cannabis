"""Preserve FNS-II alternative chemistry, protein leads and assay boundaries."""
from .phase1_reference_gap_bundle import run

SOURCES = (
    'data/curation/fnsii-carrier-review.json',
    'data/reports/phase1-fnsii-alternative-audit.json',
    'data/reports/phase1-fnsii-references.json',
    'data/reports/phase1-fnsii-search.json',
    'data/reports/phase1-fnsii-annotations.json',
    'data/curation/fnsii-rice-primary-assay-review.json',
)


if __name__ == '__main__':
    run(SOURCES, 'fnsii-evidence-bundle')
