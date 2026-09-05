"""Whole-metabolome net tests after adding thirteen screened reaction hypotheses."""
import copy
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
import scipy

from .phase1_marts_completions import balanced
from .phase1_net_flux import NetModel
from .phase1_scope import write_rows
from .phase1_screened_overlay import build_overlay


def build(network, baseline, sensitivity, plant, alternatives):
    original = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    evidence = copy.deepcopy(baseline['candidate_reaction_evidence_ids'])
    supplements = []
    for search, filename in [(plant, 'phase1-plant-purine-search.json'), (alternatives, 'phase1-purine-gap-search.json')]:
        rows = {r['reaction_id']: r for r in search['rows']}
        if len(rows) != len(search['rows']):
            raise ValueError('Duplicate search equation')
        overlay = build_overlay({'reactions': network['reactions'], 'hypotheses': []}, search, filename)
        for e in overlay['enzyme_evidence']:
            rid = e['reaction_id']
            if rid in evidence:
                raise ValueError('Supplement overlaps existing candidate equation')
            e['evidence_class'] = rows[rid]['evidence_class']
            supplements.append(e)
            evidence[rid] = [e['id']]
    if not evidence.keys() <= original.keys() or any(not ids for ids in evidence.values()):
        raise ValueError('Invalid candidate evidence inventory')
    reactions = [original[rid] for rid in sorted(evidence)]
    if any(not balanced([r['left'], r['right']], compounds) for r in reactions):
        raise ValueError('Candidate equation not element/isotope/charge balanced')
    exchange = set(baseline['external_exchange_compound_ids'])
    if exchange != set(sensitivity['external_exchange_compound_ids']):
        raise ValueError('Exchange boundary changed')
    co2 = baseline['co2_compound_id']
    if {c for c in exchange if compounds[c]['carbon_count']} != {co2} or compounds[co2]['smiles'] != 'O=C=O':
        raise ValueError('CO2 must be sole carbon exchange')
    if [(t['cannabisdb_id'], t['compound_id']) for t in baseline['targets']] != [(t['cannabisdb_id'], t['compound_id']) for t in network['targets']]:
        raise ValueError('Target inventory changed')
    old = {c['compound_id']: c for c in baseline['certificates']}
    restricted_status = {t['cannabisdb_id']: t['restricted_net_status'] for t in sensitivity['targets']}
    scenarios = []
    for name, forbidden in [('permissive-directions', set()),
                            ('five-reverse-steps-forbidden', {c['id'] for c in sensitivity['constraints']})]:
        model = NetModel(reactions, exchange, forbidden_step_ids=forbidden)
        preserved = {cid: cert for cid, cert in old.items() if not forbidden.intersection(s['step_id'] for s in cert['steps'])}
        certificates, targets, cache = {}, [], {}
        for t in baseline['targets']:
            cid = t['compound_id']
            if cid in preserved:
                result = preserved[cid]
            else:
                if cid not in cache:
                    cache[cid] = model.solve(cid)
                result = cache[cid]
            if result['status'] == 'exact-net-conversion-hypothesis':
                if not {s['reaction_id'] for s in result['steps']} <= evidence.keys() or forbidden.intersection(s['step_id'] for s in result['steps']):
                    raise ValueError('Invalid certificate step')
                carbon_in = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['external_net_consumption'].items())
                carbon_out = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['net_exports'].items())
                if carbon_in <= 0 or carbon_in != carbon_out or Fraction(result['external_net_consumption'].get(co2, '0')) != carbon_in:
                    raise ValueError('Exact CO2 carbon balance failed')
                certificates[cid] = copy.deepcopy(result) if cid in preserved else {'compound_id': cid, **result}
            prior_status = t['net_status'] if not forbidden else restricted_status[t['cannabisdb_id']]
            targets.append({k: t[k] for k in ('cannabisdb_id', 'compound_id', 'label', 'carbon_count')} | {
                'previous_scenario_status': prior_status, 'net_status': result['status'],
                'certificate_compound_id': cid if cid in certificates else None,
                'new_net_certificate': result['status'] == 'exact-net-conversion-hypothesis' and prior_status != 'exact-net-conversion-hypothesis'} |
                {k: result[k] for k in ('solver_status', 'solver_message') if k in result})
        summary = {'target_records': len(targets), 'candidate_equations': len(reactions),
            'allowed_directed_steps': len(model.steps), 'unique_reactant_compounds': len({m['compound_id'] for s in model.steps for m in s['required_inputs']}),
            'target_status_counts': dict(Counter(t['net_status'] for t in targets)),
            'exact_structure_certificates': len(certificates), 'preserved_structure_certificates': len(preserved),
            'new_target_ids': [t['cannabisdb_id'] for t in targets if t['new_net_certificate']]}
        scenarios.append({'id': name, 'summary': summary, 'forbidden_step_ids': sorted(forbidden),
            'targets': targets, 'certificates': list(certificates.values()), 'preserved_certificate_compound_ids': sorted(preserved)})
        print(name, json.dumps(summary), flush=True)
    added = {e['reaction_id'] for e in supplements}
    used = added | {s['reaction_id'] for scenario in scenarios for c in scenario['certificates'] for s in c['steps']}
    used_compounds = exchange | {m['compound_id'] for rid in used for side in ('left', 'right') for m in original[rid][side]}
    return {'schema': 'cannabis-carbon.phase1-purine-candidate-net.v1',
        'summary': {'target_records_per_scenario': len(baseline['targets']), 'candidate_equations': len(evidence),
            'added_candidate_equations': len(added), 'added_distinct_proteins': len({p['accession'] for e in supplements for p in e['screened_proteins']})},
        'scenarios': scenarios, 'candidate_reaction_evidence_ids': evidence, 'enzyme_evidence': supplements,
        'reactions': [{**original[rid], 'enzyme_evidence_ids': evidence[rid], 'new_purine_candidate': rid in added} for rid in sorted(used)],
        'compounds': [compounds[c] for c in sorted(used_compounds)],
        'external_exchange_compound_ids': sorted(exchange), 'co2_compound_id': co2,
        'constraints': copy.deepcopy(sensitivity['constraints']), 'solver': copy.deepcopy(baseline['solver']), 'scipy_version': scipy.__version__,
        'claim_boundary': 'All target records tested in two separate candidate-only scenarios; these are net conversion hypotheses with regenerated pre-existing pools, not startup or physiological pathways. All added equations have screened protein hypotheses, not confirmed activity. Reviewed and unreviewed reference evidence remain distinct. The five forbidden reverse steps are analyst-selected sensitivity constraints; all other directions remain hypothetical. No organic carbon exchange is admitted. Original reports remain unchanged. Atom tracing remains deferred.'}


def run():
    names = ('phase1-full-balanced-network', 'phase1-expanded-candidate-net', 'phase1-direction-sensitivity',
             'phase1-plant-purine-search', 'phase1-purine-gap-search')
    paths = [Path('data/reports', n + '.json') for n in names]
    reports = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for report in reports:
        for path, sha in report.get('source_sha256', {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Source lineage changed')
        if 'source_discovery' in report and hashlib.sha256(Path(report['source_discovery']).read_bytes()).hexdigest() != report['source_discovery_sha256']:
            raise ValueError('Discovery source changed')
    report = build(*reports); report['source_sha256'] = hashes
    payload = json.dumps(report, separators=(',', ':')) + '\n'; sha = hashlib.sha256(payload.encode()).hexdigest()
    Path('data/reports/phase1-purine-candidate-net.json').write_text(payload)
    groups = [('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id'), ('enzyme_evidence', 'enzyme_evidence', 'id')]
    records = [('metadata', 'report', {k: v for k, v in report.items() if k not in {'scenarios', *(g[1] for g in groups)}})]
    for kind, key, id_key in groups:
        records.extend((kind, r[id_key], r) for r in report[key])
    for scenario in report['scenarios']:
        records.append(('scenario', scenario['id'], {k: v for k, v in scenario.items() if k not in ('targets', 'certificates')}))
        for kind, key, id_key in [('target', 'targets', 'cannabisdb_id'), ('certificate', 'certificates', 'compound_id')]:
            records.extend((kind, scenario['id'] + ':' + r[id_key], r) for r in scenario[key])
    count = write_rows(records, sha, Path('data/derived/phase1-purine-candidate-net.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
