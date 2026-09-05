"""Lossless comparator and residue-review export, with no model promotion."""
from .phase1_reference_gap_bundle import run

SOURCES = (
    'data/reports/phase1-flavone-fht-comparison.json',
    'data/curation/flavone-fht-comparison-review.json',
    'data/curation/flavone-seven-site-review.json',
    'data/reports/phase1-flavone-site-review.json',
)


if __name__ == '__main__':
    run(SOURCES, 'flavone-specificity-bundle')
