"""Whole-metabolome sensitivity to exact Rhea synthase reference links."""
import copy
import hashlib
import json
from pathlib import Path

from .phase1_purine_candidate_net import build as build_net
from .phase1_scope import write_rows


def build(network, parent, links):
    scenarios = {s['id']: s for s in parent['scenarios']}
    baseline = copy.deepcopy({**parent, **scenarios['permissive-directions']})
    original_ids = set(baseline['candidate_reaction_evidence_ids'])
    constraints = copy.deepcopy(parent['constraints'])
    for link in links['rows']:
        if not link['protein_links'] or link['core_identity_merge_allowed']:
            raise ValueError('Invalid synthase evidence boundary')
        rid = link['reaction_id']
        eid = 'synthase-reference-link:' + link['id']
        baseline['candidate_reaction_evidence_ids'].setdefault(rid, []).append(eid)
        forward = link['canonical_forward_side']
        if forward not in ('left', 'right'):
            raise ValueError('Unknown reference direction')
        reverse = 'right-to-left' if forward == 'left' else 'left-to-right'
        constraints.append({'id': rid + ':hypothetical-' + reverse, 'reaction_id': rid,
            'source_reaction_id': link['annotated_rhea_id'],
            'reference_annotation': copy.deepcopy(link['reference_catalytic_annotation']),
            'reason': 'Sensitivity excludes reversal of the source-annotated synthase direction; not proof of thermodynamic irreversibility.'})
    restricted = scenarios['six-reverse-steps-forbidden']
    sensitivity = {'external_exchange_compound_ids': parent['external_exchange_compound_ids'],
        'constraints': constraints, 'targets': [{'cannabisdb_id': t['cannabisdb_id'], 'restricted_net_status': t['net_status']} for t in restricted['targets']]}
    result = build_net(network, baseline, sensitivity, None, None, search_supplements=[],
                       restricted_scenario_id='eight-reverse-steps-forbidden')
    added = set(result['candidate_reaction_evidence_ids']) - original_ids
    result['schema'] = 'cannabis-carbon.phase1-synthase-candidate-net.v1'
    result['summary']['added_candidate_equations'] = len(added)
    result['summary'].pop('added_distinct_proteins')
    result['synthase_reference_links'] = copy.deepcopy(links['rows'])
    original = {r['id']: r for r in network['reactions']}
    present = {r['id'] for r in result['reactions']}
    linked = {r['reaction_id'] for r in links['rows']}
    for rid in sorted(linked - present):
        result['reactions'].append({**copy.deepcopy(original[rid]), 'enzyme_evidence_ids': result['candidate_reaction_evidence_ids'][rid]})
    for reaction in result['reactions']:
        reaction.pop('new_purine_candidate', None)
        reaction['new_synthase_equation'] = reaction['id'] in added
    result['reactions'].sort(key=lambda r: r['id'])
    used = set(result['external_exchange_compound_ids']) | {m['compound_id'] for r in result['reactions'] for side in ('left', 'right') for m in r[side]}
    result['compounds'] = [c for c in network['compounds'] if c['id'] in used]
    result['claim_boundary'] = ('All 6,220 historical exact targets evaluated without neutralization or stereochemical merging. '
        'Adds one exact CBDAS Rhea equation and supplements existing THCAS evidence. Reference annotation and sequence identity are distinct from candidate homology and direct assays. '
        'Permissive directions remain a chemical upper bound; the restricted scenario preserves six prior reverse exclusions and adds both synthase reversals. '
        'Only CO2 supplies net carbon. Regenerated pre-existing pools are allowed; startup, expression, compartmentation and physiological flux remain unresolved. '
        'Historical identity conflicts are not repaired. Atom tracing remains deferred.')
    return result


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-full-balanced-network', 'phase1-replacement-candidate-net', 'phase1-synthase-reaction-links')]
    inputs = [json.loads(p.read_text()) for p in paths]
    for report in inputs:
        for p, sha in report.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Source lineage changed')
    report = build(*inputs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-synthase-candidate-net.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id'), ('synthase_reference_link', 'synthase_reference_links', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {'scenarios', *(g[1] for g in groups)}})]
    for kind, collection, key in groups:
        rows.extend((kind, r[key], r) for r in report[collection])
    for scenario in report['scenarios']:
        rows.append(('scenario', scenario['id'], {k: v for k, v in scenario.items() if k not in ('targets', 'certificates')}))
        for kind, collection, key in [('target', 'targets', 'cannabisdb_id'), ('certificate', 'certificates', 'compound_id')]:
            rows.extend((kind, scenario['id'] + ':' + r[key], r) for r in scenario[collection])
    count = write_rows(rows, sha, Path('data/derived/phase1-synthase-candidate-net.ndjson'))
    print(json.dumps({'summary': report['summary'], 'rows': count, 'sha256': sha}))


if __name__ == '__main__':
    run()
