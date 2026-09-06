"""Forward-only lipid acylation sensitivity on the full chemistry scenario."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_reaction_completion_net import assemble, validate_certificate
from .phase1_marts_completions import balanced
from .phase1_net_flux import NetModel


def equation_key(r):
    return tuple(sorted(tuple(sorted((m['compound_id'], m['coefficient']) for m in r[s])) for s in ('left', 'right')))


def build(network, completions, original_net, parent_net, lipid, *, prior_layers=(), extra_forward_types=()):
    reactions, compounds, _ = assemble(network, completions)
    inherited_forbidden = set()
    for layer in prior_layers:
        for c in layer['compounds']:
            if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
                raise ValueError('Prior layer compound conflict')
            compounds.setdefault(c['id'], c)
        for r in layer['added_reactions']:
            if r['id'] in reactions and equation_key(r) != equation_key(reactions[r['id']]):
                raise ValueError('Prior layer equation conflict')
            if not balanced([r['left'], r['right']], compounds):
                raise ValueError('Prior layer balance failed')
            reactions.setdefault(r['id'], r)
        inherited_forbidden.update(layer['forbidden_step_ids'])
    keys = {equation_key(r): r['id'] for r in reactions.values()}
    for c in lipid['compounds']:
        if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
            raise ValueError('Lipid compound identity conflict')
        compounds.setdefault(c['id'], c)
    added, joins, forbidden = {}, [], sorted(inherited_forbidden)
    for r in lipid['reactions']:
        if not balanced([r['left'], r['right']], compounds):
            raise ValueError('Lipid equation failed independent balance')
        key = equation_key(r)
        if key in keys:
            joins.append({'hypothesis_id': r['id'], 'existing_reaction_id': keys[key]})
            continue
        keys[key] = r['id']; reactions[r['id']] = r; added[r['id']] = r
        if r['hypothesis_type'] in ('sn1-acylation', 'sn2-acylation', 'sn3-acylation', 'phosphatidate-hydrolysis',
                                   'cardiolipin-synthesis', 'pgp-hydrolysis', 'pgp-synthesis', 'cdp-dag-synthesis',
                                   'pe-synthesis', 'ps-synthesis', 'pe-first-methylation', 'pe-second-methylation', *extra_forward_types):
            forbidden.append(r['id'] + ':hypothetical-right-to-left')
    if [(t['cannabisdb_id'], t['compound_id']) for t in network['targets']] != [(t['cannabisdb_id'], t['compound_id']) for t in parent_net['targets']]:
        raise ValueError('Target inventory mismatch')
    exchange = set(parent_net['external_exchange_compound_ids']); co2 = parent_net['co2_compound_id']
    if co2 not in exchange or compounds[co2]['smiles'] != 'O=C=O' or {c for c in exchange if compounds[c]['carbon_count']} != {co2}:
        raise ValueError('CO2 must remain the sole carbon exchange')
    if exchange != set(original_net['external_exchange_compound_ids']):
        raise ValueError('Exchange boundary changed')
    model = NetModel(list(reactions.values()), exchange, forbidden_step_ids=forbidden)
    steps = {s['id']: s for s in model.steps}
    certs = {c['compound_id']: c for c in original_net['certificates'] + parent_net['new_certificates']}
    for cert in certs.values():
        validate_certificate(cert, steps, compounds, exchange, co2)
    cache, new, targets = {}, {}, []
    participants = {m['compound_id'] for r in reactions.values() for side in ('left', 'right') for m in r[side]}
    for i, t in enumerate(parent_net['targets'], 1):
        cid = t['compound_id']
        if cid in certs:
            result = certs[cid]
        else:
            if cid not in cache:
                cache[cid] = model.solve(cid)
            result = cache[cid]
            if result['status'] == 'exact-net-conversion-hypothesis' and cid not in new:
                cert = {'compound_id': cid, **result}
                validate_certificate(cert, steps, compounds, exchange, co2)
                cert['added_lipid_reaction_ids'] = sorted({s['reaction_id'] for s in cert['steps']} & added.keys())
                if not cert['added_lipid_reaction_ids']:
                    raise ValueError('New witness lacks added lipid chemistry')
                new[cid] = cert
        targets.append({k: t[k] for k in ('cannabisdb_id', 'label', 'compound_id')} | {
            'parent_net_status': t['net_status'],
            'net_status': result['status'].replace('no-net-producing-candidate-equation', 'no-net-producing-equation'),
            'balanced_participant': cid in participants,
            'new_certificate': cid in new,
            'added_lipid_reaction_ids': new.get(cid, {}).get('added_lipid_reaction_ids', [])})
        if i % 500 == 0:
            print(f'Lipid net: {i}/{len(parent_net["targets"])}; {len(new)} new certificates', flush=True)
    used = {s['reaction_id'] for c in new.values() for s in c['steps']}
    return {'schema': 'cannabis-carbon.phase1-lipid-acylation-net.v1',
        'summary': {'target_records': len(targets), 'balanced_equations': len(reactions),
            'added_lipid_equations': len(added), 'existing_equation_joins': len(joins),
            'net_status_counts': dict(Counter(t['net_status'] for t in targets)),
            'balanced_participant_records': sum(t['balanced_participant'] for t in targets),
            'new_certificate_structures': len(new), 'new_certificate_records': sum(t['new_certificate'] for t in targets)},
        'targets': targets, 'new_certificates': list(new.values()), 'added_reactions': list(added.values()),
        'existing_equation_joins': joins, 'certificate_reactions': [reactions[r] for r in sorted(used)],
        'compounds': list(compounds.values()), 'external_exchange_compound_ids': sorted(exchange),
        'co2_compound_id': co2, 'forbidden_step_ids': sorted(forbidden),
        'claim_boundary': 'Chemistry-only sensitivity, not confirmed Cannabis pathways. New acylation equations run only in source-forward orientation; explicit PA acid-base bridges may run both ways. Existing catalog/completion equations retain their permissive directions. Exact net certificates allow regenerated pre-existing pools, not established pool origin, physiological input conditions, energetics or compartments. No enzyme or atom-mapping claim.',
        'baseline_certificate_reports': ['phase1-catalog-net-gaps.json', 'phase1-reaction-completion-net.json'],
        'lipid_evidence_report': 'phase1-lipid-acylation.json'}


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-full-balanced-network',
        'phase1-marts-completions', 'phase1-catalog-net-gaps', 'phase1-reaction-completion-net', 'phase1-lipid-acylation')]
    inputs = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for report in inputs:
        for p, digest in report.get('source_sha256', {}).items():
            if p in hashes and hashes[p] != digest:
                raise ValueError('Pinned source checksum mismatch')
    report = build(*inputs); report['source_sha256'] = hashes
    Path('data/reports/phase1-lipid-acylation-net.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
