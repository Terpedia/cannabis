"""Exact-family search for the all-target weighted diagnostic's new gaps."""
from pathlib import Path
from .phase1_weighted_gap_search import run


if __name__ == '__main__':
    run(Path('data/reports/phase1-remaining-weighted-routes.json'), 'phase1-remaining-gap')
