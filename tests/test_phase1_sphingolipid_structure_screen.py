import json
from pathlib import Path
from cannabis_carbon.phase1_sphingolipid_structure_screen import build, classify


def test_topology_screen_distinguishes_amide_and_four_hydroxy_backbone():
    assert classify('CCCC(=O)NC(CO)C(O)C(O)CCCC')['status'] == 'backbone-topology-lead'
    assert classify('CCCC(=O)NC(CO)C(O)CCCCC')['status'] == 'no-backbone-match'
    assert classify('NC(CO)C(O)C(O)CCCC')['status'] == 'no-backbone-match'


def test_entire_target_inventory_preserved(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-sphingolipid-structure-screen.json').read_text())
    rebuilt = build()
    # RDKit patch releases may differ; replay actual structural outcomes.
    rebuilt['rdkit_version'] = report['rdkit_version']
    assert report == rebuilt
    assert len(report['rows']) == 6220
    assert report['summary']['exact_structures_screened'] == 6203
    assert sum(report['summary']['target_status_counts'].values()) == 6220
    assert all(r['model_eligible'] is False for r in report['rows'])
