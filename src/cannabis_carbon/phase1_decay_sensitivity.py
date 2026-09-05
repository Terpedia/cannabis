"""Focused chemistry sensitivity, not a candidate-enzyme or completeness claim."""
import copy
import hashlib
import json
from collections import Counter
from pathlib import Path

from .phase1_net_flux import NetModel
from .phase1_scope import write_rows


def build(network, precursor, review, candidates):
    reactions = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    previous = next(s for s in precursor['scenarios'] if s['id'] == 'restricted-full-catalog-chemistry')
    exchange = set(precursor['external_exchange_compound_ids'])
    if exchange != set(candidates['external_exchange_compound_ids']):
        raise ValueError('Exchange boundary mismatch')
    co2 = precursor['co2_compound_id']
    if {c for c in exchange if compounds[c]['carbon_count']} != {co2} or compounds[co2]['smiles'] != 'O=C=O':
        raise ValueError('Organic carbon exchange forbidden')
    rid = review['reaction_id']
    join = next(j for j in reactions[rid]['sources'] if j['source_reaction_id'] == review['rhea_decay_direction'])
    if join['source_left_corresponds_to'] != 'right' or review['selected_canonical_direction'] != 'hypothetical-left-to-right':
        raise ValueError('Reviewed orientation changed')
    excluded = rid + ':' + review['selected_canonical_direction']
    constraints = copy.deepcopy(precursor['constraints'])
    if excluded in {c['id'] for c in constraints}:
        raise ValueError('Restriction already present')
    constraints.append({'id': excluded, 'reaction_id': rid, 'review_id': review['id'],
                        'boundary': 'Analyst-selected sensitivity, not demonstrated irreversibility.'})
    model = NetModel(network['reactions'], exchange, {c['id'] for c in constraints})
    results = []
    for old in previous['results']:
        cid = old['compound_id']
        result = {'compound_id': cid, **model.solve(cid), 'previous_status': old['status'],
                  'previous_used_excluded_step': any(s['step_id'] == excluded for s in old.get('steps', []))}
        results.append(result)
        print(len(results), cid, result['status'], flush=True)
    used = {s['reaction_id'] for r in results for s in r.get('steps', [])}
    needed = exchange | {r['compound_id'] for r in results} | {
        m['compound_id'] for r in used for side in ('left', 'right') for m in reactions[r][side]}
    candidate_ids = set(candidates['candidate_reaction_evidence_ids'])
    gaps = []
    for gap_id in sorted(used - candidate_ids):
        gaps.append({'reaction_id': gap_id, 'source_joins': reactions[gap_id]['sources'],
                     'selected_uses': [{'probe_compound_id': r['compound_id'], **s} for r in results
                                       for s in r.get('steps', []) if s['reaction_id'] == gap_id],
                     'status': 'no-candidate-link-in-current-1601-equation-model'})
    gaps.sort(key=lambda g: (-len(g['selected_uses']), g['reaction_id']))
    targets = {t['compound_id'] for t in precursor['focused_targets']}
    return {'schema': 'cannabis-carbon.phase1-decay-sensitivity.v1',
            'summary': {'probe_count': len(results), 'equation_count': len(reactions),
                        'allowed_directed_steps': len(model.steps),
                        'status_counts': dict(Counter(r['status'] for r in results)),
                        'focused_target_status_counts': dict(Counter(r['status'] for r in results if r['compound_id'] in targets)),
                        'previously_used_excluded_step': sum(r['previous_used_excluded_step'] for r in results),
                        'alternative_certificates': sum(r['previous_used_excluded_step'] and r['status'] == 'exact-net-conversion-hypothesis' for r in results),
                        'selected_missing_candidate_equations': len(gaps)},
            'constraints': constraints, 'excluded_step_id': excluded,
            'review': copy.deepcopy(review), 'results': results, 'candidate_gaps': gaps,
            'participant_compound_ids': precursor['participant_compound_ids'],
            'focused_targets': precursor['focused_targets'], 'co2_compound_id': co2,
            'external_exchange_compound_ids': sorted(exchange),
            'reactions': [reactions[r] for r in sorted(used)],
            'compounds': [compounds[c] for c in sorted(needed)],
            'solver': copy.deepcopy(precursor['solver']),
            'claim_boundary': 'Each of 39 exact probes is independently re-solved in the full balanced catalog with six restricted steps. Other directions remain hypothetical. Missing Cannabis protein evidence is retained, not promoted by chemical feasibility. CO2 is the sole carbon exchange; regenerated pre-existing pools remain possible. This is not joint supply, startup, thermodynamic feasibility, Cannabis activity or whole-metabolome completeness. Selected gap membership is not necessity. Historical scenarios and atom accounting are unchanged.'}


def run():
    paths = [Path(p) for p in ('data/reports/phase1-full-balanced-network.json',
             'data/reports/phase1-purine-precursor-audit.json',
             'data/curation/ureidoglycine-direction-review.json',
             'data/reports/phase1-purine-candidate-net.json')]
    inputs = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for item in inputs:
        for path, sha in item.get('source_sha256', {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Source snapshot mismatch: ' + path)
    report = build(*inputs)
    report['source_sha256'] = hashes
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-decay-sensitivity.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('result', 'results', 'compound_id'), ('reaction', 'reactions', 'id'),
              ('compound', 'compounds', 'id'), ('candidate_gap', 'candidate_gaps', 'reaction_id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, key, identifier in groups:
        rows.extend((kind, r[identifier], r) for r in report[key])
    count = write_rows(rows, sha, Path('data/derived/phase1-decay-sensitivity.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
