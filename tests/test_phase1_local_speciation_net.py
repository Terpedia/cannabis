import hashlib
import json
from collections import Counter
from pathlib import Path
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_lipid_acylation_net import equation_key
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate


def test_full_inventory_local_speciation_certificates_and_unchanged_boundaries():
    root = Path('data/reports')
    read = lambda p: json.loads(p.read_bytes())
    report = read(root/'phase1-local-speciation-net.json')
    parent = read(root/'phase1-glycerophospholipid-net.json')
    proposals = read(root/'phase1-local-speciation-hypotheses.json')
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    layers = [read(root/n) for n in report['baseline_certificate_reports']]
    reactions, compounds, model = assemble_current(read(root/'phase1-full-balanced-network.json'),
        read(root/'phase1-marts-completions.json'), report, layers)
    steps = {s['id']:s for s in model.steps}
    for k in ('forbidden_step_ids','external_exchange_compound_ids','co2_compound_id'):
        assert report[k] == parent[k]
    assert [(t['cannabisdb_id'],t['compound_id']) for t in report['targets']] == [
        (t['cannabisdb_id'],t['compound_id']) for t in parent['targets']]
    assert len(report['targets']) == report['summary']['target_records'] == 6220
    proposed = {r['id']:r for r in proposals['reactions']}
    added = {r['id']:r for r in report['added_reactions']}
    assert set(added) | {j['hypothesis_id'] for j in report['existing_equation_joins']} == set(proposed)
    for rid, r in added.items():
        assert r == proposed[rid]
        assert r['hypothesis_type'] == 'local-speciation'
        assert all(rid+':hypothetical-'+d in steps for d in ('left-to-right','right-to-left'))
    for join in report['existing_equation_joins']:
        assert equation_key(proposed[join['hypothesis_id']]) == equation_key(reactions[join['existing_reaction_id']])
    prior_certs = {c['compound_id']:c for layer in layers for c in layer.get('certificates',[]) + layer.get('new_certificates',[])}
    new = {c['compound_id']:c for c in report['new_certificates']}
    assert not new.keys() & prior_certs.keys()
    assert len(new) == len(report['new_certificates']) == report['summary']['new_certificate_structures']
    assert len(prior_certs) == 2636
    for cert in [*prior_certs.values(), *new.values()]:
        validate_certificate(cert,steps,compounds,set(report['external_exchange_compound_ids']),report['co2_compound_id'])
    for cert in new.values():
        used = {s['reaction_id'] for s in cert['steps']} & added.keys()
        assert used and used == set(cert['added_lipid_reaction_ids'])
    parent_targets = {t['cannabisdb_id']:t for t in parent['targets']}
    participants = {p['compound_id'] for r in reactions.values() for side in ('left','right') for p in r[side]}
    for target in report['targets']:
        assert target['parent_net_status'] == parent_targets[target['cannabisdb_id']]['net_status']
        assert target['new_certificate'] == (target['compound_id'] in new)
        assert target['balanced_participant'] == (target['compound_id'] in participants)
        if target['compound_id'] in prior_certs or target['compound_id'] in new:
            assert target['net_status'] == 'exact-net-conversion-hypothesis'
    assert report['summary']['net_status_counts'] == dict(Counter(t['net_status'] for t in report['targets']))
    assert report['summary']['new_certificate_records'] == sum(t['new_certificate'] for t in report['targets'])
    assert report['summary']['balanced_participant_records'] == sum(t['balanced_participant'] for t in report['targets'])
    assert report['summary']['balanced_equations'] == parent['summary']['balanced_equations'] + len(added)
