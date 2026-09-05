import json
from pathlib import Path

from cannabis_carbon.phase1_weak_hit_domain_review import build


def test_domain_review_preserves_sequences_and_rejects_catalytic_promotion(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-weak-hit-domain-review.json').read_text())
    assert report == build()
    assert report['decision']['model_eligible'] is False
    assert report['coordinates'] == {'qstart': 8, 'qend': 89, 'sstart': 18, 'send': 99}
    for protein in report['proteins']:
        assert any(x['id'] == 'PF00173' for x in protein['domain_crossrefs'])
        domain = next(f for f in protein['features'] if f['type'] == 'Domain')
        assert domain['description'] == 'Cytochrome b5 heme-binding'
    assert report['alignment']['passes_screen'] is True
