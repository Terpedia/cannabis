"""Publish cofactor-specific gaps, broad reference searches, and domain review."""
from .phase1_reference_gap_bundle import run

SOURCES = tuple('data/reports/phase1-' + name + '.json' for name in (
    'ureide-plant-reference-review', 'ureide-nonplant-reference-review',
    'ureidoglycolate-broad-references', 'ureidoglycolate-broad-search',
    'ureidoglycolate-domain-review')) + ('data/curation/ureidoglycolate-cofactor-review.json',)


if __name__ == '__main__':
    run(SOURCES, 'ureidoglycolate-evidence-bundle')
