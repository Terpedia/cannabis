"""Separate chemical sensitivity excluding two unsupported ureide condensations."""
from .phase1_allantoate_sensitivity import run as solve

RID = 'balanced-equation:2365a3fda8ba5acf89d328e4e34a798ee5c268e3c00261cd77506ec26376419d'


def run():
    solve(additional_forbidden=[RID + ':hypothetical-left-to-right'],
          name='ureide-sensitivity',
          review_path='data/curation/ureide-dual-condensation-review.json')


if __name__ == '__main__':
    run()
