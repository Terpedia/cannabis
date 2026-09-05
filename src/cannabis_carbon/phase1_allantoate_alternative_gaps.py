"""Identify exact missing-enzyme steps in condensation-free witnesses."""
import hashlib
import json
from pathlib import Path


def build(search_paths=None, route_path=Path('data/reports/phase1-allantoate-sensitivity.json')):
    route_path = Path(route_path)
    model_path = Path('data/reports/phase1-remaining-candidate-net.json')
    route, model = [json.loads(p.read_text()) for p in (route_path, model_path)]
    paths = [route_path, model_path]
    prior = {}
    for path in search_paths if search_paths is not None else sorted(Path('data/reports').glob('phase1-*search.json')):
        path = Path(path)
        report = json.loads(path.read_text())
        if not isinstance(report.get('rows'), list):
            continue
        paths.append(path)
        for row in report['rows']:
            if 'reaction_id' in row and 'search_status' in row:
                prior.setdefault(row['reaction_id'], []).append({'report': str(path), 'row': row})
    reactions = {r['id']: r for r in route['reactions']}
    gaps = {}
    for target in route['rows']:
        for step in target.get('steps', []):
            rid = step['reaction_id']
            if rid in model['candidate_reaction_evidence_ids']:
                continue
            gap = gaps.setdefault(rid, {'reaction_id': rid, 'target_ids': [], 'selected_uses': [],
                'source_joins': reactions[rid]['sources'], 'prior_searches': prior.get(rid, [])})
            gap['target_ids'].append(target['cannabisdb_id'])
            gap['selected_uses'].append({**step, 'target_ids': [target['cannabisdb_id']], 'compound_id': target['compound_id']})
    for gap in gaps.values():
        gap['target_ids'] = sorted(set(gap['target_ids']))
    rows = sorted(gaps.values(), key=lambda r: (-len(r['target_ids']), r['reaction_id']))
    return {'schema': 'cannabis-allantoate-alternative-gaps-v1', 'candidate_gaps': rows,
        'reactions': route['reactions'], 'compounds': route['compounds'],
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'missing_enzyme_equations': len(rows), 'previously_unsearched_equations': sum(not r['prior_searches'] for r in rows)},
        'claim_boundary': 'Missing-enzyme steps in fourteen condensation-free chemical witnesses, not supported biosynthesis or guaranteed rescues. Search history is pinned to the listed report snapshots, including reviewed rejections; presence of a historical search is not enzyme evidence.'}


if __name__ == '__main__':
    report = build()
    Path('data/reports/phase1-allantoate-alternative-gaps.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))
