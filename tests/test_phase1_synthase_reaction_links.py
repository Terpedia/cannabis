import copy
import hashlib
import json
from pathlib import Path

from cannabis_carbon.phase1_synthase_reaction_links import build
from cannabis_carbon.balance import _reaction_smiles_balance

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_exact_rhea_joins_preserve_annotation_direction_and_all_sequence_links():
    inputs = [read(n) for n in ('phase1-full-balanced-network', 'phase1-cannabinoid-revalidation-references', 'phase1-cannabinoid-revalidation-search')]
    before = copy.deepcopy(inputs)
    report = read('phase1-synthase-reaction-links')
    assert build(*inputs) == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert inputs == before
    assert report['summary'] == {'exact_rhea_equations': 2, 'protein_reaction_links': 67, 'sequence_identical_links': 1, 'historical_core_equations_merged': 0}
    assert {r['annotated_rhea_id'] for r in report['rows']} == {'RHEA:34412', 'RHEA:34136'}
    for row in report['rows']:
        assert row['core_identity_merge_allowed'] is False
        assert row['reaction_id'] != row['core_reaction_id']
        assert row['canonical_forward_side'] == 'left'
        sides = {c['role']: c for c in row['core_comparison']}
        assert not sides['substrates']['exact_encoded_match']
        assert not sides['substrates']['uncharger_only_match']
        assert sides['substrates']['uncharger_and_stereo_removed_match']
        assert not sides['products']['exact_encoded_match']
        assert sides['products']['uncharger_only_match']
        expected = {a['id'] for a in inputs[2]['passing_alignments'] if a['reference_accession'] == row['reference_accession']}
        assert {p['alignment_id'] for p in row['protein_links']} == expected
        assert all(not p['direct_candidate_assay_claimed'] for p in row['protein_links'])
    exact = [p for r in report['rows'] for p in r['protein_links'] if p['sequence_identical']]
    assert exact[0]['candidate_accession'] == 'A0A7J6G9C8'
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_all_four_distinct_equations_are_balanced_without_identity_merging():
    report = read('phase1-synthase-reaction-links')
    compounds = {c['id']: c for c in report['compounds']}
    assert len(report['reactions']) == 4
    for reaction in report['reactions']:
        sides = ['.'.join(compounds[m['compound_id']]['smiles'] for m in reaction[side] for _ in range(m['coefficient'])) for side in ('left', 'right')]
        element, charge = _reaction_smiles_balance('>>'.join(sides))
        assert element['status'] == charge['status'] == 'balanced'
