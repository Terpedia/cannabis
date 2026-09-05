"""Publish desaturase rejections and alternative-substrate hypotheses losslessly."""
from .phase1_reference_gap_bundle import run

SOURCES = tuple('data/reports/phase1-' + name + '.json' for name in (
    'weak-plant-reference-review', 'weak-nonplant-reference-review', 'weak-nonplant-search',
    'weak-hit-domain-review', 'reviewed-weak-search', 'long-desaturase-review',
    'desaturase-domain-references', 'desaturase-domain-search',
    'sphingolipid-alternative-references', 'sphingolipid-alternative-search',
    'sphingolipid-catalog-review', 'sphingolipid-structure-screen'))

if __name__ == '__main__':
    run(SOURCES, 'desaturase-evidence-bundle')
