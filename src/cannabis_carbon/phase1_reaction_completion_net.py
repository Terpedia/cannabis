"""Chemistry-only sensitivity: full catalog plus balanced completion hypotheses."""
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced
from .phase1_net_flux import NetModel, exact_net
from .phase1_scope import expand

BOUNDARY = ('Reaction-completeness sensitivity, not confirmed Cannabis pathways. '
    'No enzyme gate. Added equations are inferred stoichiometric completions, not '
    'validated source corrections. Both directions are hypothetical. CO2 is the '
    'only carbon exchange; the unchanged carbon-free reservoir is permissive. '
    'Net conversion may require pre-existing regenerated pools; startup, '
    'physiological direction, thermodynamics and compartments remain unresolved. '
    'Atom tracing is deferred. Baseline evidence is not promoted.')


def assemble(network, completions):
    compounds = {c['id']: c for c in network['compounds']}
    reactions = {r['id']: r for r in network['reactions']}
    added = {}
    for c in completions['compounds']:
        if stable_id('structure', c['smiles']) != c['id']:
            raise ValueError('Compound identity digest mismatch')
        if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
            raise ValueError('Compound identity conflict')
        compounds.setdefault(c['id'], c)
    for h in completions['completions']:
        rid = h['balanced_equation_id']
        if stable_id('balanced-equation', [h['left'], h['right']]) != rid:
            raise ValueError('Equation identity mismatch')
        if rid in reactions:
            if any(h[s] != reactions[rid][s] for s in ('left', 'right')):
                raise ValueError('Equation join conflict')
            continue
        added.setdefault(rid, {'id': rid, 'left': h['left'], 'right': h['right'],
            'completion_ids': [], 'evidence_class': 'inferred-stoichiometric-completion',
            'direction_status': 'both-directions-hypothetical', 'enzyme_evidence_ids': []})['completion_ids'].append(h['id'])
    reactions.update(added)
    for r in reactions.values():
        if not balanced([r['left'], r['right']], compounds):
            raise ValueError('Equation failed element/isotope/charge balance')
    return reactions, compounds, added


def validate_certificate(cert, steps, compounds, exchange, co2):
    selected = [steps[s['step_id']] for s in cert['steps']]
    for s, original in zip(cert['steps'], selected):
        if s['reaction_id'] != original['reaction_id'] or s['direction_mode'] != original['direction_mode']:
            raise ValueError('Certificate step identity mismatch')
    net = exact_net(selected, [s['extent'] for s in cert['steps']])
    if net.get(cert['compound_id'], 0) < 1 or any(n < 0 for c, n in net.items() if c not in exchange):
        raise ValueError('Certificate depletes an internal input or fails target production')
    incoming = {c: str(-n) for c, n in sorted(net.items()) if n < 0}
    outgoing = {c: str(n) for c, n in sorted(net.items()) if n > 0}
    if incoming != cert['external_net_consumption'] or outgoing != cert['net_exports']:
        raise ValueError('Certificate net accounting mismatch')
    carbon_in = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in incoming.items())
    carbon_out = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in outgoing.items())
    if carbon_in <= 0 or carbon_in != carbon_out or Fraction(incoming.get(co2, '0')) != carbon_in:
        raise ValueError('CO2 carbon balance failed')


def build(network, completions, baseline):
    reactions, compounds, added = assemble(network, completions)
    if [(t['cannabisdb_id'], t['compound_id']) for t in network['targets']] != [(t['cannabisdb_id'], t['compound_id']) for t in baseline['targets']]:
        raise ValueError('Target inventory mismatch')
    exchange = set(baseline['external_exchange_compound_ids'])
    co2 = baseline['co2_compound_id']
    if co2 not in exchange or compounds[co2]['smiles'] != 'O=C=O' or {c for c in exchange if compounds[c]['carbon_count']} != {co2}:
        raise ValueError('CO2 must be the sole carbon exchange')
    model = NetModel(list(reactions.values()), exchange)
    steps = {s['id']: s for s in model.steps}
    startup = expand(list(reactions.values()), exchange)
    certs = {c['compound_id']: c for c in baseline['certificates']}
    for c in certs.values():
        validate_certificate(c, steps, compounds, exchange, co2)
    cache, targets, new = {}, [], {}
    participants = {m['compound_id'] for r in reactions.values() for side in ('left', 'right') for m in r[side]}
    for i, target in enumerate(baseline['targets'], 1):
        cid = target['compound_id']
        if cid in certs:
            result = certs[cid]
        else:
            if cid not in cache:
                cache[cid] = model.solve(cid)
            result = cache[cid]
            if result['status'] == 'exact-net-conversion-hypothesis' and cid not in new:
                cert = {'compound_id': cid, **result}
                validate_certificate(cert, steps, compounds, exchange, co2)
                cert['added_reaction_ids'] = sorted({s['reaction_id'] for s in cert['steps']} & added.keys())
                if not cert['added_reaction_ids']:
                    raise ValueError('New certificate has no added chemistry; review baseline comparison')
                new[cid] = cert
        status = result['status'].replace('no-net-producing-candidate-equation', 'no-net-producing-equation')
        targets.append({k: target[k] for k in ('cannabisdb_id', 'label', 'compound_id')} | {
            'baseline_net_status': target['net_status'], 'net_status': status,
            'balanced_participant': cid in participants,
            'has_net_producer_in_hypothetical_direction': cid in model.producible,
            'startup_status': 'external-exchange' if cid in exchange else 'structural-scope-reachable' if cid in startup['available'] else 'blocked',
            'new_certificate': cid in new,
            'added_reaction_ids': new.get(cid, {}).get('added_reaction_ids', [])})
        if i % 500 == 0:
            print(f'Reaction completion: {i}/{len(baseline["targets"])}; {len(new)} new exact certificates', flush=True)
    used = {s['reaction_id'] for c in new.values() for s in c['steps']}
    return {'schema': 'cannabis-carbon.phase1-reaction-completion-net.v1', 'claim_boundary': BOUNDARY,
        'summary': {'target_records': len(targets), 'baseline_equations': len(network['reactions']),
            'added_completion_equations': len(added), 'total_balanced_equations': len(reactions),
            'net_status_counts': dict(Counter(t['net_status'] for t in targets)),
            'balanced_participant_records': sum(t['balanced_participant'] for t in targets),
            'records_without_net_producer': sum(not t['has_net_producer_in_hypothetical_direction'] for t in targets),
            'new_certificate_structures': len(new), 'new_certificate_records': sum(t['new_certificate'] for t in targets),
            'startup_status_counts': dict(Counter(t['startup_status'] for t in targets))},
        'targets': targets, 'new_certificates': list(new.values()),
        'added_reactions': list(added.values()), 'certificate_reactions': [reactions[r] for r in sorted(used)],
        'compounds': list(compounds.values()), 'external_exchange_compound_ids': sorted(exchange),
        'co2_compound_id': co2, 'baseline_certificate_report': 'phase1-catalog-net-gaps.json',
        'completion_evidence_report': 'phase1-marts-completions.json'}


def run():
    paths = [Path('data/reports', n + '.json') for n in
        ('phase1-full-balanced-network', 'phase1-marts-completions', 'phase1-catalog-net-gaps')]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    inputs = [json.loads(p.read_text()) for p in paths]
    for report in inputs:
        for p, digest in report.get('source_sha256', {}).items():
            if p in hashes and hashes[p] != digest:
                raise ValueError('Pinned source checksum mismatch')
    report = build(*inputs)
    report['source_sha256'] = hashes
    Path('data/reports/phase1-reaction-completion-net.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    export_report()
    print(json.dumps(report['summary']), flush=True)


def export_report(name='reaction-completion-net'):
    path = Path('data/reports/phase1-' + name + '.json')
    payload = path.read_bytes()
    report = json.loads(payload)
    digest = hashlib.sha256(payload).hexdigest()
    Path('docs/data/' + name + '.json').write_bytes(payload)
    rows = [{'record_kind': 'source_document', 'record_id': str(path),
             'record_json': json.dumps(report, separators=(',', ':')), 'report_sha256': digest}]
    Path('data/derived/phase1-' + name + '.ndjson').write_text(
        ''.join(json.dumps(r, separators=(',', ':')) + '\n' for r in rows))


if __name__ == '__main__':
    run()
