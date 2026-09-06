import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from cannabis_carbon.phase1_c17_evidence_audit import annotate
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate


def test_annotation_separates_source_presence_direction_and_hypotheses():
    compounds = {'a': {'smiles': 'OO'}, 'b': {'smiles': 'O=O'}}
    step = {'required_inputs': [{'compound_id': 'a'}], 'outputs': [{'compound_id': 'b'}],
            'direction_mode': 'hypothetical-left-to-right'}
    r = {'id': 'r', 'sources': [{'source_urls': ['https://example.invalid/source']}], 'hypothesis_type': 'analogy'}
    a = annotate(r, step, compounds)
    assert a['evidence_class'] == 'proposed-gap-filling-reaction'
    assert a['review_flags'] == ['peroxide-consuming-oxygen-producing-direction-review']
    assert not a['cannabis_physiological_direction_established_by_this_audit']
    assert not a['enzyme_assignment_established_by_this_audit']
    assert annotate({'id': 'r'}, step, compounds)['evidence_class'] == 'source-link-or-stoichiometric-completion-review-required'
    reverse = {**step, 'required_inputs': step['outputs'], 'outputs': step['required_inputs']}
    assert annotate(r, reverse, compounds)['review_flags'] == []


def test_all_five_saved_certificates_evidence_and_inputs_replay():
    root = Path('data/reports')
    read = lambda n: json.loads((root / ('phase1-' + n + '.json')).read_bytes())
    report = read('c17-evidence-audit'); current = read('c17-elongation-net')
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    layers = [json.loads((root / n).read_bytes()) for n in current['baseline_certificate_reports']]
    reactions, compounds, model = assemble_current(read('full-balanced-network'), read('marts-completions'), current, layers)
    steps = {s['id']: s for s in model.steps}
    precursor = next(t for t in read('odd-chain-supply')['targets'] if t['acyl_carbon_count'] == 15)
    expected = {c['compound_id']: c for c in current['new_certificates']}
    expected[precursor['compound_id']] = {'compound_id': precursor['compound_id'], **precursor['net_result']}
    assert len(report['certificates']) == len(expected) == 5
    assert {c['cannabisdb_id'] for c in report['certificates']} == {None, 'CDB000153', 'CDB000155', 'CDB000157', 'CDB000386'}
    used = set()
    for row in report['certificates']:
        cert = row['certificate']
        assert cert == expected[row['compound_id']]
        validate_certificate(cert, steps, compounds, set(current['external_exchange_compound_ids']), current['co2_compound_id'])
        assert len(row['steps']) == len(cert['steps'])
        for item, saved in zip(row['steps'], cert['steps']):
            step = steps[saved['step_id']]
            assert item == {**saved, **annotate(reactions[saved['reaction_id']], step, compounds),
                            'required_inputs': step['required_inputs'], 'outputs': step['outputs']}
            used.update(p['compound_id'] for side in ('required_inputs', 'outputs') for p in step[side])
        assert row['summary'] == {'step_count': len(row['steps']),
            'evidence_class_counts': dict(Counter(s['evidence_class'] for s in row['steps'])),
            'review_flag_step_counts': dict(Counter(f for s in row['steps'] for f in s['review_flags']))}
        assert {i['compound_id']: Fraction(i['amount_per_target']) for i in row['net_inputs']} == {
            c: Fraction(v) / Fraction(cert['target_amount']) for c, v in cert['external_net_consumption'].items()}
        for i in row['net_inputs']:
            assert i['smiles'] == compounds[i['compound_id']]['smiles']
            assert not i['nutrient_uptake_established']
    assert report['compounds'] == [compounds[c] for c in sorted(used)]
    assert report['summary'] == {'audited_certificates': 5, 'new_inventory_target_certificates': 4,
        'new_enzyme_assignments': 0, 'new_confirmed_Cannabis_pathways': 0, 'coverage_change': 0}
