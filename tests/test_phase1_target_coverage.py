from cannabis_carbon.phase1_target_coverage import encoded_structure, audit_targets
import hashlib
import json
from pathlib import Path


def test_exact_encoding_does_not_merge_stereo_charge_or_salts():
    assert encoded_structure('C[C@H](O)Cl')[0] != encoded_structure('C[C@@H](O)Cl')[0]
    assert encoded_structure('C[O-]')[0] != encoded_structure('CO')[0]
    assert encoded_structure('CO.[Na+]')[0] != encoded_structure('CO')[0]
    assert encoded_structure('[13CH4]')[0] != encoded_structure('C')[0]
    assert encoded_structure('CC(O)CC=CC')[1] == 'stereo-unspecified-or-unknown'


def test_every_target_retained_and_xref_does_not_grant_structural_coverage():
    targets = [{'id': 'A', 'smiles': 'CCO'}, {'id': 'B', 'smiles': None},
               {'id': 'C', 'smiles': 'CC', 'external_ids': {'chebi': '1'}}, {'id': 'D', 'smiles': 'CO'}]
    network = {'reactions': [
        {'id': 'rhea:1', 'reaction_smiles': 'C=C.O>>CCO', 'reactants': [], 'products': [{'compound_id': 'chebi:1'}]},
        {'id': 'rhea:2', 'reaction_smiles': 'CO.*>>CC', 'reactants': [], 'products': []}]}
    result = audit_targets(targets, network, {'reactions': []})
    assert len(result['targets']) == 4
    a, b, c, d = result['targets']
    assert a['coverage_status'] == 'balanced-reaction-participant'
    assert a['balanced_right_side_record_count'] == 1
    assert b['coverage_status'] == 'structure-unresolved'
    assert c['coverage_status'] == 'reaction-participant-balance-unresolved'
    assert c['xref_reaction_records'] == ['core:rhea:1']
    assert c['balanced_reaction_record_count'] == 0
    assert d['coverage_status'] == 'reaction-participant-balance-unresolved'


def test_published_target_audit_keeps_entire_inventory_and_valid_reaction_links():
    root = Path(__file__).resolve().parents[1]
    source = root / 'docs/data/compounds.json'
    targets = json.loads(source.read_text())['compounds']
    report = json.loads((root / 'docs/data/phase1-target-coverage.json').read_text())
    assert report['source_sha256']['docs/data/compounds.json'] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert [x['id'] for x in targets] == [x['cannabisdb_id'] for x in report['targets']]
    reactions = {r['id']: r for r in report['reaction_ledger']}
    for row in report['targets']:
        balanced = [x for x in row['reaction_matches'] if x['computed_balance_status'] == 'balanced']
        assert row['balanced_reaction_record_count'] == len(balanced)
        assert (row['coverage_status'] == 'balanced-reaction-participant') == bool(balanced)
        for match in row['reaction_matches']:
            record = reactions[match['reaction_record_id']]
            assert record['computed_balance_status'] == match['computed_balance_status']
            for role in match['roles']:
                side = record['reaction_smiles'].split('>>')[0 if role['equation_side'] == 'left' else 1]
                matches = sum(encoded_structure(fragment)[0] == row['canonical_isomeric_smiles'] for fragment in side.split('.'))
                assert matches == role['coefficient']
