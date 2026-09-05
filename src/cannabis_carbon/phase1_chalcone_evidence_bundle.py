"""Lossless CHI reference, search, annotation and residue-review export."""
from .phase1_reference_gap_bundle import run

SOURCES = ('data/curation/chalcone-reference-review.json',) + tuple(
    'data/reports/phase1-' + name + '.json' for name in (
        'chalcone-references', 'chalcone-search', 'chalcone-annotations', 'chalcone-site-review'))


if __name__ == '__main__':
    run(SOURCES, 'chalcone-evidence-bundle')
