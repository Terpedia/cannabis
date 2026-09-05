"""Re-solve candidate-only CO2 routes after catalog/backfill protein discovery.

This is distinct from annotating a frozen chemistry-only certificate. Every
admitted equation needs an existing candidate evidence link, but such a link
does not establish activity, physiological direction or complex assembly.
"""
import copy
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import scipy

from .phase1_marts_completions import balanced
from .phase1_net_flux import NetModel
from .phase1_scope import expand, write_rows


def build(network, startup, supplement, baseline):
    compounds = {c['id']: c for c in network['compounds']}
    original = {r['id']: r for r in network['reactions']}
    evidence = copy.deepcopy(startup['candidate_reaction_evidence_ids'])
    added = {}
    for e in supplement['enzyme_evidence']:
        rid = e['reaction_id']
        if rid in added or rid in evidence or rid not in original:
            raise ValueError('Duplicate, baseline or unknown added equation')
        if not e['screened_proteins']:
            raise ValueError('Added equation lacks screened proteins')
        added[rid] = e['id']
        evidence[rid] = [e['id']]
    if not evidence.keys() <= original.keys() or any(not ids for ids in evidence.values()):
        raise ValueError('Candidate evidence inventory mismatch')
    reactions = [{**copy.deepcopy(original[rid]), 'enzyme_evidence_ids': ids,
                  'new_catalog_candidate': rid in added} for rid, ids in sorted(evidence.items())]
    for r in reactions:
        if not balanced([r['left'], r['right']], compounds):
            raise ValueError('Candidate equation fails element/isotope/charge balance')
    exchanges = set(baseline['external_exchange_compound_ids'])
    co2 = baseline['co2_compound_id']
    carbon_exchanges = {cid for cid in exchanges if compounds[cid]['carbon_count']}
    if carbon_exchanges != {co2} or compounds[co2]['smiles'] != 'O=C=O':
        raise ValueError('CO2 must be the sole carbon exchange')
    old_scenario = next(s for s in startup['scenarios'] if s['id'] == 'CO2-plus-all-carbon-free-species')
    if set(old_scenario['seed_compound_ids']) != exchanges:
        raise ValueError('Startup exchange boundary changed')
    if [t['cannabisdb_id'] for t in baseline['targets']] != [t['cannabisdb_id'] for t in network['targets']]:
        raise ValueError('Target inventory changed')
    scope = expand(reactions, exchanges)
    if not set(old_scenario['witnesses']) <= scope['available']:
        raise ValueError('Additive scope lost baseline compounds')
    model = NetModel(reactions, exchanges)
    old_certificates = {c['compound_id']: c for c in baseline['certificates']}
    certificates, targets, cache = {}, [], {}
    for i, target in enumerate(baseline['targets'], 1):
        cid = target['compound_id']
        if cid in old_certificates:
            # Addition cannot invalidate a valid witness. Preserve it verbatim.
            result = old_certificates[cid]
        else:
            if cid not in cache:
                cache[cid] = model.solve(cid)
            result = cache[cid]
        if result['status'] == 'exact-net-conversion-hypothesis':
            carbon_in = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['external_net_consumption'].items())
            carbon_out = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['net_exports'].items())
            if carbon_in <= 0 or carbon_in != carbon_out or Fraction(result['external_net_consumption'].get(co2, '0')) != carbon_in:
                raise ValueError('Exact CO2 carbon balance failed')
            if not {s['reaction_id'] for s in result['steps']} <= evidence.keys():
                raise ValueError('Certificate uses an unsupported equation')
            certificates[cid] = copy.deepcopy(result) if cid in old_certificates else {
                'compound_id': cid, **result, 'net_carbon_in': str(carbon_in), 'net_carbon_out': str(carbon_out)}
        targets.append({k: target[k] for k in ('cannabisdb_id', 'compound_id', 'label', 'carbon_count')} | {
            'baseline_net_status': target['net_status'], 'net_status': result['status'],
            'baseline_startup_status': target['startup_status'],
            'startup_status': 'explicit-seed' if cid in exchanges else 'hypothesis-scope-reachable' if cid in scope['available'] else 'blocked',
            'certificate_compound_id': cid if cid in certificates else None,
            'new_net_certificate': cid in certificates and cid not in old_certificates} |
            {k: result[k] for k in ('solver_status', 'solver_message') if k in result})
        if i % 500 == 0:
            print(f'Expanded candidates: {i}/{len(baseline["targets"])} targets; {len(certificates)} certificates', flush=True)
    used = {s['reaction_id'] for c in certificates.values() for s in c['steps']}
    used_compounds = {m['compound_id'] for r in reactions if r['id'] in used for side in ('left', 'right') for m in r[side]} | exchanges
    new_ids = [t['cannabisdb_id'] for t in targets if t['new_net_certificate']]
    return {'schema': 'cannabis-carbon.phase1-expanded-candidate-net.v1',
        'summary': {'target_records': len(targets), 'candidate_equations': len(reactions),
            'added_candidate_equations': len(added), 'exact_structure_certificates': len(certificates),
            'new_structure_certificates': len(certificates.keys() - old_certificates.keys()),
            'new_target_certificates': len(new_ids), 'new_target_ids': new_ids,
            'target_status_counts': dict(Counter(t['net_status'] for t in targets)),
            'startup_status_counts': dict(Counter(t['startup_status'] for t in targets)),
            'new_startup_target_ids': [t['cannabisdb_id'] for t in targets if t['baseline_startup_status'] == 'blocked' and t['startup_status'] == 'hypothesis-scope-reachable'],
            'selected_new_candidate_equations': len(used & added.keys())},
        'targets': targets, 'certificates': list(certificates.values()),
        'reactions': [r for r in reactions if r['id'] in used],
        'compounds': [compounds[c] for c in sorted(used_compounds)],
        'candidate_reaction_evidence_ids': evidence, 'enzyme_evidence': supplement['enzyme_evidence'],
        'startup_witnesses': scope['witnesses'],
        'external_exchange_compound_ids': sorted(exchanges), 'co2_compound_id': co2,
        'solver': baseline['solver'], 'scipy_version': scipy.__version__,
        'claim_boundary': baseline['claim_boundary'] + ' Candidate-only expansion, not the frozen catalog certificate annotation count. All previous successful certificates are preserved verbatim; other targets are re-solved with the 116 added exact-equation candidate links. No inferred MARTS completions or archived-reference extensions are included. One homolog does not establish required complex assembly or full catalytic capability.'}


def run():
    paths = [Path('data/reports', name + '.json') for name in ('phase1-full-balanced-network',
        'phase1-candidate-scope', 'phase1-combined-catalog-evidence', 'phase1-candidate-net-flux')]
    reports = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for report in reports:
        for path, sha in report.get('source_sha256', {}).items():
            if path in hashes and hashes[path] != sha:
                raise ValueError('Expanded candidate source lineage mismatch')
    output = build(*reports)
    output['source_sha256'] = hashes
    payload = json.dumps(output, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-expanded-candidate-net.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('target', 'targets', 'cannabisdb_id'), ('certificate', 'certificates', 'compound_id'),
              ('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id'), ('enzyme_evidence', 'enzyme_evidence', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in output.items() if k not in {g[1] for g in groups}})]
    for kind, collection, key in groups:
        rows.extend((kind, r[key], r) for r in output[collection])
    count = write_rows(rows, sha, Path('data/derived/phase1-expanded-candidate-net.ndjson'))
    print(json.dumps({'summary': output['summary'], 'rows': count, 'sha256': sha}), flush=True)


if __name__ == '__main__':
    run()
