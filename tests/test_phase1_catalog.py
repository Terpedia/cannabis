import pytest
import json
from pathlib import Path
from cannabis_carbon.phase1_catalog import assemble
from cannabis_carbon.balance import _reaction_smiles_balance


def fixture():
    base = {'source_urls': ['https://example.org/source'], 'source_uniprot_ids': [], 'source_genbank_ids': [], 'source_ec_numbers': []}
    good = {**base, 'reaction_id': 'RHEA:1', 'reaction_smarts': 'C=C.O>>CCO', 'balance_status': 'balanced'}
    bad = {**base, 'reaction_id': 'MARTS:1', 'reaction_smarts': 'C=C>>CCO', 'balance_status': 'imbalanced', 'source_uniprot_ids': ['P1']}
    alt = {'reaction_id': bad['reaction_id'], 'original_reaction_smarts': bad['reaction_smarts'],
           'balanced_reference_candidates': [{'rule_id': good['reaction_id'], 'reaction_smarts': good['reaction_smarts'],
             'carbon_stoichiometry_matches': True, 'participant_changes': [{'smiles': 'O', 'coefficient_delta': 1}], 'required_review': 'verify source'}]}
    return {'rows': [good, bad]}, {'rows': [good]}, {'rows': [alt]}, {'rows': []}


def test_alternatives_are_links_not_extra_reactions_or_enzyme_transfers():
    result = assemble(*fixture())
    assert result['summary']['balanced_reaction_variants'] == 1
    assert len(result['source_ledger']) == 2
    assert len(result['compounds']) == 3  # Includes carbon-free water.
    reaction = result['reactions'][0]
    assert len(reaction['reactants']) == 2
    assert reaction['source_references']['source_uniprot_ids'] == []
    assert reaction['alternative_source_links'][0]['original_source_references']['source_uniprot_ids'] == ['P1']
    assert result['source_ledger'][1]['catalog_reaction_id'] is None


def test_missing_ledger_or_false_balance_fails_closed():
    queue, evidence, alternatives, family = fixture()
    with pytest.raises(ValueError, match='every excluded'):
        assemble(queue, evidence, {'rows': []}, family)
    queue['rows'][0]['reaction_smarts'] = 'C=C>>CCO'
    with pytest.raises(ValueError, match='fails the current audit'):
        assemble(queue, evidence, alternatives, family)


def test_full_catalog_preserves_participants_coefficients_and_all_original_records():
    root = Path(__file__).resolve().parents[1]
    report = json.loads((root / 'docs/data/phase1-reaction-catalog.json').read_text())
    queue = json.loads((root / 'data/reports/phase1-enzyme-discovery-queue.json').read_text())
    compounds = {row['id']: row for row in report['compounds']}
    reactions = {row['id']: row for row in report['reactions']}
    assert len(compounds) == len(report['compounds'])
    assert len(reactions) == len(report['reactions'])
    assert {(r['reaction_id'], r['reaction_smarts']) for r in queue['rows']} == {
        (r['reaction_id'], r['reaction_smarts']) for r in report['source_ledger']}
    for row in report['source_ledger']:
        if row['catalog_reaction_id']:
            assert row['original_balance_status'] == 'balanced'
            assert row['catalog_reaction_id'] in reactions
        for alternative in row['balanced_alternative_links']:
            assert alternative['catalog_reaction_id'] in reactions
    for reaction in reactions.values():
        sides = []
        for side in ('reactants', 'products'):
            fragments = []
            for participant in reaction[side]:
                assert isinstance(participant['coefficient'], int) and participant['coefficient'] > 0
                fragments.extend([compounds[participant['compound_id']]['smiles']] * participant['coefficient'])
            sides.append('.'.join(fragments))
        element, charge = _reaction_smiles_balance('>>'.join(sides))
        assert element['differences'] == {}
        assert charge['difference'] == 0
