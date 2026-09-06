import hashlib
import json
from collections import Counter
from pathlib import Path
import pytest
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_net_obstruction import validate as validate_obstruction


@pytest.mark.parametrize('report_name,parent_name,proposal_name,added_count,equations', [
    ('odd-chain-net','alkane-net','odd-chain-elongation',8,18211),
    ('c17-elongation-net','odd-chain-net','c17-elongation',4,18215)])
def test_full_inventory_exact_replay_and_changed_obstructions(report_name,parent_name,proposal_name,added_count,equations):
    root = Path('data/reports')
    read = lambda n: json.loads((root / ('phase1-' + n + '.json')).read_bytes())
    report = read(report_name); parent = read(parent_name)
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    layers = [json.loads((root / n).read_bytes()) for n in report['baseline_certificate_reports']]
    _, compounds, model = assemble_current(read('full-balanced-network'), read('marts-completions'), report, layers)
    steps = {s['id']: s for s in model.steps}
    prior = {c['compound_id']: c for layer in layers for c in layer.get('certificates', []) + layer.get('new_certificates', [])}
    new = {c['compound_id']: c for c in report['new_certificates']}
    assert len(prior) == 2729
    assert len(new) == len(report['new_certificates']) and not prior.keys() & new.keys()
    for c in [*prior.values(), *new.values()]:
        validate_certificate(c, steps, compounds, set(parent['external_exchange_compound_ids']), parent['co2_compound_id'])
    added = {r['id']: r for r in report['added_reactions']}
    assert added == {r['id']: r for r in read(proposal_name)['reactions']}
    assert len(added) == added_count and not report['existing_equation_joins']
    for rid in added:
        assert rid + ':hypothetical-left-to-right' in steps
        assert rid + ':hypothetical-right-to-left' not in steps
    assert set(report['forbidden_step_ids']) == set(parent['forbidden_step_ids']) | {rid + ':hypothetical-right-to-left' for rid in added}
    assert report['external_exchange_compound_ids'] == parent['external_exchange_compound_ids']
    assert [(t['cannabisdb_id'], t['compound_id']) for t in report['targets']] == [(t['cannabisdb_id'], t['compound_id']) for t in parent['targets']]
    assert len(report['targets']) == 6220
    assert report['summary']['balanced_equations'] == equations
    assert report['summary']['net_status_counts'] == dict(Counter(t['net_status'] for t in report['targets']))
    for t in report['targets']:
        assert (t['net_status'] == 'exact-net-conversion-hypothesis') == (t['compound_id'] in prior or t['compound_id'] in new)
        assert t['new_certificate'] == (t['compound_id'] in new)
    for c in new.values():
        used = {s['reaction_id'] for s in c['steps']} & added.keys()
        assert used and set(c['added_lipid_reaction_ids']) == used
    if report_name == 'c17-elongation-net':
        assert {t['cannabisdb_id'] for t in report['targets'] if t['new_certificate']} == {
            'CDB000153', 'CDB000155', 'CDB000157', 'CDB000386'}
        assert len(new) == 4
        assert all(set(c['added_lipid_reaction_ids']) == set(added) for c in new.values())
    # The prior proofs belong to the prior model, not this expanded one.
    for obstruction in read('odd-chain-obstructions' if report_name == 'c17-elongation-net' else 'alkane-obstructions')['targets']:
        with pytest.raises(ValueError, match='Allowed step increases obstruction weight'):
            validate_obstruction(model, obstruction['compound_id'], obstruction['weights'])
