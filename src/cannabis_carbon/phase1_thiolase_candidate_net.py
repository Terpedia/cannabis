"""Whole-target and exact synthase probes after the weighted-gap protein search."""
import copy
import hashlib
import json
from pathlib import Path

from .phase1_net_flux import NetModel
from .phase1_purine_candidate_net import build as build_net
from .phase1_scope import write_rows


def build(network, parent, search, links, *, search_filename='phase1-weighted-gap-search.json'):
    scenarios = {s['id']: s for s in parent['scenarios']}
    baseline = {**parent, **scenarios['permissive-directions']}
    restricted_id = 'eight-reverse-steps-forbidden'
    sensitivity = {'external_exchange_compound_ids': parent['external_exchange_compound_ids'],
        'constraints': parent['constraints'], 'targets': [
            {'cannabisdb_id': t['cannabisdb_id'], 'restricted_net_status': t['net_status']}
            for t in scenarios[restricted_id]['targets']]}
    result = build_net(network, baseline, sensitivity, None, None,
        search_supplements=[(search, search_filename)],
        restricted_scenario_id=restricted_id)
    result['schema'] = 'cannabis-carbon.phase1-thiolase-candidate-net.v1'
    result['synthase_reference_links'] = copy.deepcopy(parent['synthase_reference_links'])
    result['probes'] = copy.deepcopy(links['compounds'])
    selected = [r for r in network['reactions'] if r['id'] in result['candidate_reaction_evidence_ids']]
    probe_results = []
    for scenario in result['scenarios']:
        model = NetModel(selected, result['external_exchange_compound_ids'], scenario['forbidden_step_ids'])
        for compound in result['probes']:
            cid = compound['id']
            row = {'id': scenario['id'] + ':' + cid, 'scenario_id': scenario['id'],
                'compound_id': cid, **model.solve(cid)}
            probe_results.append(row)
            print('probe', scenario['id'], cid, row['status'], flush=True)
    result['probe_results'] = probe_results
    used = {r['id'] for r in result['reactions']} | {s['reaction_id'] for r in probe_results for s in r.get('steps', [])}
    added = set(result['candidate_reaction_evidence_ids']) - parent['candidate_reaction_evidence_ids'].keys()
    result['reactions'] = [{**copy.deepcopy(r), 'enzyme_evidence_ids': result['candidate_reaction_evidence_ids'][r['id']],
        'new_thiolase_candidate': r['id'] in added} for r in network['reactions'] if r['id'] in used]
    needed = set(result['external_exchange_compound_ids']) | {c['id'] for c in result['probes']} | {
        m['compound_id'] for r in result['reactions'] for side in ('left', 'right') for m in r[side]}
    result['compounds'] = [c for c in network['compounds'] if c['id'] in needed]
    result['claim_boundary'] = ('All 6,220 historical exact targets retained and tested separately from eight exact synthase participant probes. '
        'The added equation has reviewed-reference homology candidates, not characterized Cannabis activity. '
        'Eight previous reverse exclusions are preserved; all other orientations remain hypotheses, including acetyl-CoA condensation. '
        'CO2 is the only net carbon input. Net certificates allow regenerated pre-existing pools, not zero-pool startup. '
        'No neutralization or stereochemical merging; compartmentation, expression, physiological flux and atom tracing remain unresolved.')
    return result


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-full-balanced-network',
        'phase1-synthase-candidate-net', 'phase1-weighted-gap-search', 'phase1-synthase-reaction-links')]
    inputs = [json.loads(p.read_text()) for p in paths]
    for report in inputs:
        for p, sha in report.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Source lineage changed')
        if 'source_discovery' in report and hashlib.sha256(Path(report['source_discovery']).read_bytes()).hexdigest() != report['source_discovery_sha256']:
            raise ValueError('Discovery lineage changed')
    report = build(*inputs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-thiolase-candidate-net.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id'),
        ('enzyme_evidence', 'enzyme_evidence', 'id'), ('synthase_reference_link', 'synthase_reference_links', 'id'),
        ('probe_result', 'probe_results', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {'scenarios', *(g[1] for g in groups)}})]
    for kind, collection, key in groups:
        rows.extend((kind, r[key], r) for r in report[collection])
    for scenario in report['scenarios']:
        rows.append(('scenario', scenario['id'], {k: v for k, v in scenario.items() if k not in ('targets', 'certificates')}))
        for kind, collection, key in [('target', 'targets', 'cannabisdb_id'), ('certificate', 'certificates', 'compound_id')]:
            rows.extend((kind, scenario['id'] + ':' + r[key], r) for r in scenario[collection])
    count = write_rows(rows, sha, Path('data/derived/phase1-thiolase-candidate-net.ndjson'))
    print(json.dumps({'summary': report['summary'], 'rows': count, 'sha256': sha}), flush=True)


if __name__ == '__main__':
    run()
