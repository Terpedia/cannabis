import hashlib
import json
from collections import Counter
from pathlib import Path
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate


def test_full_inventory_geranial_reduction_certificates():
    root = Path('data/reports')
    read = lambda n: json.loads((root / ('phase1-' + n + '.json')).read_bytes())
    report = read('geranial-reduction-net'); parent = read('alcohol-acetates-net')
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    layers = [json.loads((root / n).read_bytes()) for n in report['baseline_certificate_reports']]
    reactions, compounds, model = assemble_current(read('full-balanced-network'), read('marts-completions'), report, layers)
    steps = {s['id']: s for s in model.steps}
    prior = {c['compound_id']: c for layer in layers for c in layer.get('certificates', []) + layer.get('new_certificates', [])}
    new = {c['compound_id']: c for c in report['new_certificates']}
    assert len(prior) == 2735 and len(new) == len(report['new_certificates']) == 2
    assert not prior.keys() & new.keys()
    exchange = set(parent['external_exchange_compound_ids']); co2 = parent['co2_compound_id']
    assert {c for c in exchange if compounds[c]['carbon_count']} == {co2}
    for c in [*prior.values(), *new.values()]:
        validate_certificate(c, steps, compounds, exchange, co2)
    added = {r['id']: r for r in report['added_reactions']}
    assert added == {r['id']: r for r in read('geranial-reduction')['reactions']}
    assert len(added) == 1 and not report['existing_equation_joins']
    for rid in added:
        assert rid + ':hypothetical-left-to-right' in steps
        assert rid + ':hypothetical-right-to-left' not in steps
    assert set(report['forbidden_step_ids']) == set(parent['forbidden_step_ids']) | {rid + ':hypothetical-right-to-left' for rid in added}
    assert report['external_exchange_compound_ids'] == parent['external_exchange_compound_ids']
    assert [(t['cannabisdb_id'], t['compound_id']) for t in report['targets']] == [(t['cannabisdb_id'], t['compound_id']) for t in parent['targets']]
    assert len(report['targets']) == 6220 and len(reactions) == report['summary']['balanced_equations'] == 18222
    assert report['summary']['net_status_counts'] == dict(Counter(t['net_status'] for t in report['targets']))
    assert report['summary']['net_status_counts']['exact-net-conversion-hypothesis'] == 2740
    assert {t['cannabisdb_id'] for t in report['targets'] if t['new_certificate']} == {'CDB000081', 'CDB000585'}
    participants = {p['compound_id'] for r in reactions.values() for side in ('left', 'right') for p in r[side]}
    for t in report['targets']:
        assert (t['net_status'] == 'exact-net-conversion-hypothesis') == (t['compound_id'] in prior or t['compound_id'] in new)
        assert t['new_certificate'] == (t['compound_id'] in new)
        assert t['balanced_participant'] == (t['compound_id'] in participants)
    for c in new.values():
        used = {s['reaction_id'] for s in c['steps']} & added.keys()
        assert used and set(c['added_lipid_reaction_ids']) == used
    assert {r['id']: r for r in report['certificate_reactions']} == {rid: reactions[rid] for c in new.values() for rid in {s['reaction_id'] for s in c['steps']}}
