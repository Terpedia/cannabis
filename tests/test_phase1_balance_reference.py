from cannabis_carbon.phase1_balance_reference import carbon_participants, match_references
import hashlib
import json
from pathlib import Path
from cannabis_carbon.balance import _reaction_smiles_balance


def test_exact_identity_keeps_stereo_charge_and_isotopes():
    assert carbon_participants('C[C@H](O)Cl>>C') != carbon_participants('C[C@@H](O)Cl>>C')
    assert carbon_participants('C[O-]>>C') != carbon_participants('CO>>C')
    assert carbon_participants('[13CH4]>>C') != carbon_participants('C>>C')
    assert carbon_participants('[CH4:1]>>C') == carbon_participants('C>>C')
    assert carbon_participants('*C>>C') is None


def test_source_equation_retains_missing_water_and_coefficient_changes():
    audit = {'reactions': [
        {'reaction_id': 'gap1', 'reaction_smarts': 'C=C>>CCO', 'status': 'imbalanced'},
        {'reaction_id': 'gap2', 'reaction_smarts': 'C=C>>CCCC', 'status': 'imbalanced'},
        {'reaction_id': 'gap3', 'reaction_smarts': '*>>C', 'status': 'not_auditable'}]}
    catalog = [
        {'rule_id': 'RHEA:1', 'reaction_smarts': 'C=C.O>>CCO'},
        {'rule_id': 'RHEA:2', 'reaction_smarts': 'C=C>>CCO'},
        {'rule_id': 'RHEA:3', 'reaction_smarts': 'C=C.C=C.[H][H]>>CCCC'}]
    rows = match_references(audit, catalog)
    assert len(rows) == 3
    first = rows[0]['balanced_reference_candidates']
    assert len(first) == 1
    assert first[0]['carbon_stoichiometry_matches'] is True
    assert '.O>>' in first[0]['reaction_smarts']
    assert first[0]['participant_changes'] == [{'side': 'reactants', 'coefficient_delta': 1, 'smiles': 'O', 'formula': 'H2O', 'formal_charge': 0, 'carbon_count': 0}]
    assert rows[1]['balanced_reference_candidates'][0]['carbon_stoichiometry_matches'] is False
    assert any(p['coefficient_delta'] == 1 and p['smiles'] == 'C=C' for p in rows[1]['balanced_reference_candidates'][0]['participant_changes'])
    assert rows[2]['balanced_reference_candidates'] == []


def test_published_alternatives_are_balanced_and_preserve_every_original_gap():
    root = Path(__file__).resolve().parents[1]
    audit_path = root / 'data/reports/terpene-identity-set-candidate-expansion-balance-audit.json'
    audit = json.loads(audit_path.read_text())
    report = json.loads((root / 'docs/data/phase1-balance-reference.json').read_text())
    assert report['source_audit_sha256'] == hashlib.sha256(audit_path.read_bytes()).hexdigest()
    assert {(r['reaction_id'], r['original_reaction_smarts']) for r in report['rows']} == {
        (r['reaction_id'], r['reaction_smarts']) for r in audit['reactions'] if r['status'] != 'balanced'}
    for row in report['rows']:
        original = carbon_participants(row['original_reaction_smarts'])
        for candidate in row['balanced_reference_candidates']:
            element, charge = _reaction_smiles_balance(candidate['reaction_smarts'])
            assert element['differences'] == {}
            assert charge['difference'] == 0
            assert [set(side) for side in original] == [set(side) for side in candidate['carbon_participants']]
            assert candidate['carbon_stoichiometry_matches'] == (original == candidate['carbon_participants'])
            assert candidate['participant_changes']
