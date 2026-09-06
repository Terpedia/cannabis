"""Lossless exact-redox derivation and factorial CO2-route sensitivity."""
from .phase1_reference_gap_bundle import run

SOURCES = (
    'data/reports/phase1-fnsii-redox-hypothesis.json',
    'data/reports/phase1-fnsii-addition-sensitivity.json',
)


if __name__ == '__main__':
    run(SOURCES, 'fnsii-route-sensitivity-bundle')
