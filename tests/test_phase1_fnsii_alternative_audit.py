import hashlib
import json
from collections import Counter
from pathlib import Path
import pytest
from cannabis_carbon.phase1_fnsii_alternative_audit import IDS, inspect_components


def test_generic_components_are_not_silently_concrete():
    result = inspect_components('*C.O>>*C=O')
    assert result['left'][0]['dummy_atoms'] == 1
    assert result['left'][0]['explicit_carbon_atoms'] == 1
    assert len(result['left']) == 2
    with pytest.raises(ValueError):
        inspect_components('C')


def test_alternative_preserves_catalog_exclusions_carriers_and_parent_gap(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-fnsii-alternative-audit.json').read_text())
    network = json.loads(Path('data/reports/phase1-full-balanced-network.json').read_text())
    catalog = json.loads(Path('data/raw/phase1-balance-reference-catalog.json').read_text())
    parent = json.loads(Path('data/reports/phase1-chalcone-remaining-gaps.json').read_text())
    assert report['parent_fnsi_gap'] == parent['candidate_gaps'][0]
    assert report['network_exclusions'] == [r for r in network['excluded_rhea_source_records'] if r['source_reaction_id'] in IDS]
    assert not any(s['source_reaction_id'] in IDS for r in network['reactions'] for s in r['sources'])
    assert [r['record'] for r in report['source_records']] == [r for r in catalog if r['rule_id'] in IDS]
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    by_id = {r['record']['rule_id']: r for r in report['source_records']}
    for row in by_id.values():
        assert row['component_audit'] == inspect_components(row['record']['reaction_smarts'])
        assert all(sum(c['dummy_atoms'] for c in row['component_audit'][side]) == 10 for side in ('left', 'right'))
    forward = by_id['RHEA:57681']['component_audit']
    reverse = by_id['RHEA:57682']['component_audit']
    for a, b in [('left', 'right'), ('right', 'left')]:
        assert Counter(c['source_smiles'] for c in forward[a]) == Counter(c['source_smiles'] for c in reverse[b])
    review = report['review']
    assert review == json.loads(Path('data/curation/fnsii-carrier-review.json').read_text())
    assert review['reduced_carrier']['participant'] == 'RHEA-COMP:11964'
    assert review['oxidized_carrier']['participant'] == 'RHEA-COMP:11965'
    assert {r['accession'] for r in report['reference_leads']} == {'Q0JFI2', 'E9KBR8', 'Q9XGT9'}
    assert all(r['model_eligible'] is False for r in report['reference_leads'])
    assert report['model_eligible'] is False
    assert report['summary']['new_exact_reactions'] == report['summary']['new_exact_enzyme_assignments'] == 0
