"""Lossless publication of ureide direction sensitivities and negative searches."""
from .phase1_reference_gap_bundle import run

SOURCES = tuple('data/reports/phase1-' + name + '.json' for name in (
    'current-gap-priority', 'allantoate-sensitivity', 'allantoate-alternative-gaps',
    'ureide-sensitivity', 'ureide-alternative-gaps', 'ureide-gap-references',
    'ureide-gap-search')) + ('data/curation/ureide-dual-condensation-review.json',)


if __name__ == '__main__':
    run(SOURCES, 'ureide-evidence-bundle')
