"""Whole-metabolome test of two new enzyme hypotheses with citrate sensitivity."""
import copy
import hashlib
import json
from pathlib import Path

from .phase1_purine_candidate_net import build as build_net
from .phase1_scope import write_rows


def build(network, parent, search, discovery):
    scenarios = {s['id']: s for s in parent['scenarios']}
    baseline = {**parent, **scenarios['permissive-directions']}
    restricted = scenarios['five-reverse-steps-forbidden']
    source_rows = {r['reaction_id']: r for r in discovery['rows']}
    if len(source_rows) != len(discovery['rows']) or len(search['rows']) != len(source_rows):
        raise ValueError('Incomplete or duplicate discovery inventory')
    for row in search['rows']:
        if row['reference_matches'] != source_rows[row['reaction_id']]['reference_matches']:
            raise ValueError('Reference annotation join changed')
    citrate = next(r for r in discovery['rows'] if 'RHEA:16846' in r['source_reaction_ids'])
    if citrate['rhea_families']['RHEA:16846']['RHEA_ID_LR'] != 'RHEA:16846':
        raise ValueError('Published citrate direction family changed')
    source = next(s for s in citrate['sources'] if s['source_reaction_id'] == 'RHEA:16846')
    if source['source_left_corresponds_to'] != 'right':
        raise ValueError('Citrate canonical orientation changed')
    constraint = {'id': citrate['reaction_id'] + ':hypothetical-left-to-right',
                  'reaction_id': citrate['reaction_id'], 'source_reaction_id': 'RHEA:16846',
                  'source_joins': copy.deepcopy(citrate['sources']),
                  'reason': 'Exclude reverse citrate-synthase usage as a sensitivity; homology does not establish citrate cleavage.',
                  'claim_boundary': 'Analyst-selected restriction, not proven physiological irreversibility.'}
    sensitivity = {'external_exchange_compound_ids': parent['external_exchange_compound_ids'],
                   'constraints': [*copy.deepcopy(parent['constraints']), constraint],
                   'targets': [{'cannabisdb_id': t['cannabisdb_id'], 'restricted_net_status': t['net_status']} for t in restricted['targets']]}
    report = build_net(network, baseline, sensitivity, None, None,
                       search_supplements=[(search, 'phase1-replacement-search.json')],
                       restricted_scenario_id='six-reverse-steps-forbidden')
    report['schema'] = 'cannabis-carbon.phase1-replacement-candidate-net.v1'
    for reaction in report['reactions']:
        reaction['new_replacement_candidate'] = reaction.pop('new_purine_candidate')
    report['claim_boundary'] = ('All 6,220 CannabisDB target records evaluated in separate candidate-only net scenarios. '
        'Two sequence-supported equations are added to the prior 1,601-equation model; neither is characterized Cannabis activity. '
        'The restricted scenario preserves five prior constraints and additionally excludes reverse citrate-synthase use. '
        'Ureidoglycine decay remains absent from this candidate model, not silently promoted as spontaneous. '
        'Other directions remain hypothetical. Preserved exact certificates are retained; every other target is reconsidered with the augmented model. '
        'CO2 is the only carbon exchange. Regenerated pre-existing pools, startup, compartmentation, thermodynamics and physiological flux remain unresolved. '
        'Prior evidence, including unreviewed-reference hypotheses and withheld partial-reference leads, is not reclassified. Historical scenarios remain unchanged. Atom tracing is deferred.')
    return report


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-full-balanced-network',
             'phase1-purine-candidate-net', 'phase1-replacement-search', 'phase1-replacement-references')]
    reports = [json.loads(p.read_text()) for p in paths]
    for report in reports:
        for path, sha in report.get('source_sha256', {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Source lineage changed')
    if reports[2]['source_discovery_sha256'] != hashlib.sha256(paths[3].read_bytes()).hexdigest():
        raise ValueError('Discovery source changed')
    report = build(*reports)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-replacement-candidate-net.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id'), ('enzyme_evidence', 'enzyme_evidence', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {'scenarios', *(g[1] for g in groups)}})]
    for kind, key, identifier in groups:
        rows.extend((kind, r[identifier], r) for r in report[key])
    for scenario in report['scenarios']:
        rows.append(('scenario', scenario['id'], {k: v for k, v in scenario.items() if k not in ('targets', 'certificates')}))
        for kind, key, identifier in [('target', 'targets', 'cannabisdb_id'), ('certificate', 'certificates', 'compound_id')]:
            rows.extend((kind, scenario['id'] + ':' + r[identifier], r) for r in scenario[key])
    count = write_rows(rows, sha, Path('data/derived/phase1-replacement-candidate-net.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
