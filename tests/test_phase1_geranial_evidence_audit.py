import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate


def test_four_gains_full_sources_direction_flags_and_exact_inputs():
    root = Path('data/reports')
    read = lambda n: json.loads((root / ('phase1-' + n + '.json')).read_bytes())
    audit = read('geranial-evidence-audit'); current = read('geranial-reduction-net')
    for p, sha in audit['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    layers = [json.loads((root / n).read_bytes()) for n in current['baseline_certificate_reports']]
    reactions, compounds, model = assemble_current(read('full-balanced-network'), read('marts-completions'), current, layers)
    steps = {s['id']: s for s in model.steps}
    expected = {c['compound_id']: c for c in current['new_certificates'] + read('alcohol-acetates-net')['new_certificates']}
    assert len(audit['certificates']) == len(expected) == 4
    assert {c['cannabisdb_id'] for c in audit['certificates']} == {'CDB000081', 'CDB000110', 'CDB000585', 'CDB000695'}
    used = set(); external_supported = 0
    for row in audit['certificates']:
        cert = row['certificate']; assert cert == expected[row['compound_id']]
        validate_certificate(cert, steps, compounds, set(current['external_exchange_compound_ids']), current['co2_compound_id'])
        assert len(row['steps']) == len(cert['steps'])
        for item, saved in zip(row['steps'], cert['steps']):
            step = steps[saved['step_id']]; reaction = reactions[saved['reaction_id']]
            for key, value in saved.items():
                assert item[key] == value
            assert item['source_record'] == reaction
            assert item['source_evidence_type'] == reaction.get('source_evidence_type')
            for side in ('required_inputs', 'outputs'):
                assert item[side] == step[side]
                used.update(p['compound_id'] for p in step[side])
            ins = {compounds[p['compound_id']]['smiles'] for p in step['required_inputs']}
            outs = {compounds[p['compound_id']]['smiles'] for p in step['outputs']}
            flags = []
            if 'OO' in ins and 'O=O' in outs:
                flags.append('peroxide-consuming-oxygen-producing-direction-review')
            if 'O=C=O' in ins:
                flags.append('co2-consuming-step-not-automatically-photosynthetic-fixation')
            if reaction.get('hypothesis_type') == 'geranial-S-citronellal':
                flags.append('external-organism-evidence-not-Cannabis-enzyme-assignment')
                assert item['evidence_class'] == 'substrate-specific-cascade-biochemistry-outside-Cannabis'
                assert item['biochemical_evidence'] == read('geranial-reduction')['biochemical_evidence']
                external_supported += 1
            elif reaction.get('hypothesis_type'):
                assert item['evidence_class'] == 'proposed-gap-filling-reaction'
            elif reaction.get('sources'):
                assert item['evidence_class'] == 'catalog-linked-equation-not-Cannabis-activity-proof'
            else:
                assert item['evidence_class'] == 'source-link-or-stoichiometric-completion-review-required'
            assert item['review_flags'] == flags
            assert not item['cannabis_physiological_direction_established_by_this_audit']
            assert not item['enzyme_assignment_established_by_this_audit']
        assert row['summary'] == {'step_count': len(row['steps']),
            'evidence_class_counts': dict(Counter(s['evidence_class'] for s in row['steps'])),
            'review_flag_step_counts': dict(Counter(f for s in row['steps'] for f in s['review_flags']))}
        assert {i['compound_id']: Fraction(i['amount_per_target']) for i in row['net_inputs']} == {
            c: Fraction(v) / Fraction(cert['target_amount']) for c, v in cert['external_net_consumption'].items()}
        for i in row['net_inputs']:
            assert i['smiles'] == compounds[i['compound_id']]['smiles']
            assert i['carbon_count'] == compounds[i['compound_id']]['carbon_count']
            assert not i['nutrient_uptake_established']
    assert external_supported == 2
    assert audit['compounds'] == [compounds[c] for c in sorted(used)]
    assert audit['summary'] == {'audited_certificates': 4, 'new_enzyme_assignments': 0,
        'new_confirmed_Cannabis_pathways': 0, 'coverage_change': 0}
