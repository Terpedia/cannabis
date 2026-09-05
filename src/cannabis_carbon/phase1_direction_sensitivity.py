"""Disable the five reviewed reverse steps, retaining every other hypothesis."""
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
import scipy
from .phase1_net_flux import NetModel
from .phase1_scope import write_rows


def build(network, expanded, review):
    constraints = []
    for r in review['reviews']:
        source_left = r['source_left_corresponds_to']
        if source_left not in ('left', 'right'):
            raise ValueError('Unknown source correspondence')
        reverse = 'hypothetical-right-to-left' if source_left == 'left' else 'hypothetical-left-to-right'
        constraints.append({'id': r['reaction_id'] + ':' + reverse, 'reaction_id': r['reaction_id'],
            'forbidden_direction_mode': reverse, 'review_id': r['id'], 'source_master_id': r['source_master_id'],
            'basis': 'analyst-selected sensitivity restriction, not a physiological direction assertion'})
    forbidden = {c['id'] for c in constraints}
    if len(forbidden) != len(constraints):
        raise ValueError('Duplicate direction restriction')
    evidence = expanded['candidate_reaction_evidence_ids']
    reactions = [r for r in network['reactions'] if r['id'] in evidence]
    if {r['id'] for r in reactions} != evidence.keys():
        raise ValueError('Candidate inventory mismatch')
    exchange = set(expanded['external_exchange_compound_ids'])
    compounds = {c['id']: c for c in network['compounds']}
    if {c for c in exchange if compounds[c]['carbon_count']} != {expanded['co2_compound_id']}:
        raise ValueError('Carbon exchange changed')
    model = NetModel(reactions, exchange, forbidden_step_ids=forbidden)
    old = {c['compound_id']: c for c in expanded['certificates']}
    preserved = {cid: c for cid, c in old.items() if not forbidden & {s['step_id'] for s in c['steps']}}
    targets, cache, alternatives = [], {}, {}
    for target in expanded['targets']:
        cid = target['compound_id']
        if cid in preserved:
            result = preserved[cid]
        else:
            if cid not in cache:
                cache[cid] = model.solve(cid)
            result = cache[cid]
        if result['status'] == 'exact-net-conversion-hypothesis' and cid not in preserved:
            if forbidden & {s['step_id'] for s in result['steps']}:
                raise ValueError('Forbidden step in alternative certificate')
            carbon_in = sum(Fraction(n) * compounds[c]['carbon_count'] for c,n in result['external_net_consumption'].items())
            carbon_out = sum(Fraction(n) * compounds[c]['carbon_count'] for c,n in result['net_exports'].items())
            if carbon_in <= 0 or carbon_in != carbon_out or Fraction(result['external_net_consumption'].get(expanded['co2_compound_id'],'0')) != carbon_in:
                raise ValueError('Alternative carbon balance failed')
            alternatives[cid] = {'compound_id': cid, **result}
        targets.append({k: target[k] for k in ('cannabisdb_id','compound_id','label','carbon_count')} | {
            'unrestricted_net_status': target['net_status'], 'restricted_net_status': result['status'],
            'was_new_expanded_target': target['new_net_certificate'],
            'certificate_origin': 'unchanged-expanded-witness' if cid in preserved else 'new-restricted-witness' if cid in alternatives else None,
            'certificate_compound_id': cid if cid in preserved or cid in alternatives else None} |
            {k:result[k] for k in ('solver_status','solver_message') if k in result})
    return {'schema':'cannabis-carbon.phase1-direction-sensitivity.v1',
        'constraints':constraints, 'targets':targets, 'alternative_certificates':list(alternatives.values()),
        'preserved_certificate_compound_ids':sorted(preserved),
        'external_exchange_compound_ids':sorted(exchange), 'co2_compound_id':expanded['co2_compound_id'],
        'summary':{'target_records':len(targets), 'candidate_equations':len(reactions),
            'forbidden_directed_steps':len(forbidden), 'allowed_directed_steps':len(model.steps),
            'unique_allowed_reactant_compounds':len({p['compound_id'] for s in model.steps for p in s['required_inputs']}),
            'target_status_counts':dict(Counter(t['restricted_net_status'] for t in targets)),
            'new_alternative_structure_certificates':len(alternatives),
            'new_expanded_targets_retaining_certificates':[t['cannabisdb_id'] for t in targets if t['was_new_expanded_target'] and t['certificate_compound_id']]},
        'solver':expanded['solver'], 'scipy_version':scipy.__version__,
        'claim_boundary':'Joint five-direction sensitivity restriction only, not a curated physiological network. Other directions, carbon-free exchanges and regenerated pre-existing pools remain permissive. Every target is retained. Existing witnesses are reused only after verifying they contain no forbidden step; all other target structures are re-solved. Solver-reported infeasibility is not proof of biological absence or mathematical infeasibility. No atom tracing. The unrestricted 108-certificate scenario is preserved.'}


def run():
    paths = [Path('data/reports', n+'.json') for n in ('phase1-full-balanced-network','phase1-expanded-candidate-net','phase1-candidate-direction-review')]
    reports = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for r in reports:
        for path,sha in r.get('source_sha256',{}).items():
            if path in hashes and hashes[path] != sha:
                raise ValueError('Sensitivity lineage mismatch')
    report = build(*reports); report['source_sha256'] = hashes
    payload = json.dumps(report,separators=(',',':'))+'\n'
    Path('data/reports/phase1-direction-sensitivity.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('constraint','constraints','id'),('target','targets','cannabisdb_id'),('certificate','alternative_certificates','compound_id')]
    rows = [('metadata','report',{k:v for k,v in report.items() if k not in {g[1] for g in groups}})]
    for kind,key,id_key in groups:
        rows.extend((kind,r[id_key],r) for r in report[key])
    count = write_rows(rows,sha,Path('data/derived/phase1-direction-sensitivity.ndjson'))
    print(json.dumps({'summary':report['summary'],'rows':count,'sha256':sha}),flush=True)


if __name__ == '__main__':
    run()
