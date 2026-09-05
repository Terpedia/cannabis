"""Lossless publication of DOPA-lyase experimental leads and site caveats."""
from .phase1_reference_gap_bundle import run

SOURCES = ('data/curation/dopa-lyase-reference-review.json',) + tuple(
    'data/reports/phase1-' + name + '.json' for name in (
        'dopa-lyase-references', 'dopa-lyase-search',
        'dopa-candidate-annotations', 'dopa-site-review'))


if __name__ == '__main__':
    run(SOURCES, 'dopa-evidence-bundle')
