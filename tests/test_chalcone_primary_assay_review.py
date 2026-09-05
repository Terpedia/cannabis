import json
from pathlib import Path
import pytest


def test_reference_assay_keeps_units_substrate_and_non_cannabis_boundary(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    review = json.loads(Path('data/curation/chalcone-primary-assay-review.json').read_text())
    old = json.loads(Path('data/curation/chalcone-reference-review.json').read_text())
    assert review['reaction_id'] == old['reaction_id']
    assert review['model_eligible'] is False
    assert review['assayed_organism'] == 'Medicago sativa'
    assay = review['assay']
    assert assay['substrate_chebi'] == 'CHEBI:15413'
    assert assay['spontaneous_background_corrected'] is True
    assert assay['kcat_per_min']['value'] / assay['k_uncatalyzed_per_min'] == pytest.approx(assay['reported_rate_enhancement'], rel=0.01)
    protein = json.loads(Path('data/raw/chalcone-annotations/P28012.json').read_text())
    constants = [m for c in protein['comments'] for m in c.get('kineticParameters', {}).get('michaelisConstants', [])]
    exact = [m for m in constants if m['substrate'].startswith(assay['substrate'])]
    assert len(exact) == 1
    assert exact[0]['constant'] == assay['Km_uM']['value']
    assert any(e.get('id') == review['source']['pmid'] for e in exact[0]['evidences'])
