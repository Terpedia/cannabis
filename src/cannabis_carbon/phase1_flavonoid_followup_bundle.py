"""Preserve genomic checks, assay review and conditional flavonoid witnesses."""
from .phase1_reference_gap_bundle import run

SOURCES = tuple('data/reports/phase1-' + name + '.json' for name in (
    'chalcone-gene-model', 'chalcone-genomic-translation', 'chalcone-addition-sensitivity',
    'chalcone-remaining-gaps', 'flavone-references', 'flavone-search', 'flavone-annotations')) + (
    'data/curation/chalcone-primary-assay-review.json', 'data/curation/flavone-specificity-review.json')


if __name__ == '__main__':
    run(SOURCES, 'flavonoid-followup-bundle')
