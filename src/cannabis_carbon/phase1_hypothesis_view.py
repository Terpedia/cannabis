"""Static, lazy-loaded Cytoscape bundles derived from the verified hypothesis catalog."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def generate(report_path, output):
    raw = report_path.read_bytes()
    report = json.loads(raw)
    compounds = {c['id']: c for c in report['compounds']}
    evidence = {e['id']: e for e in report['enzyme_evidence']}
    targets = {t['cannabisdb_id']: t for t in report['targets']}
    labels = defaultdict(list)
    for target in targets.values():
        if target['compound_id']:
            labels[target['compound_id']].append(target['label'])
    by_reaction, by_target, shards = defaultdict(list), defaultdict(list), defaultdict(dict)
    reactions = {r['id']: r for r in report['reactions']}
    for hypothesis in report['hypotheses']:
        by_reaction[hypothesis['reaction_id']].append(hypothesis)
        by_target[hypothesis['cannabisdb_id']].append({
            'id': hypothesis['id'], 'reaction_id': hypothesis['reaction_id'],
            'source_reaction_ids': sorted({s['source_reaction_id'] for s in reactions[hypothesis['reaction_id']]['sources']}),
            'has_candidate_enzyme_evidence': hypothesis['has_candidate_enzyme_evidence'],
            'net_target_coefficient': hypothesis['net_target_coefficient']})
    for rid, hypotheses in by_reaction.items():
        reaction = reactions[rid]
        cids = {m['compound_id'] for side in ('left', 'right') for m in reaction[side]}
        eids = {eid for h in hypotheses for eid in h['evidence_ids']}
        shards[rid.split(':')[1][:2]][rid] = {
            'reaction': reaction, 'hypotheses': hypotheses,
            'compounds': [{**compounds[cid], 'labels': labels[cid]} for cid in sorted(cids)],
            'enzyme_evidence': [evidence[eid] for eid in sorted(eids)]}
    def write(relative, payload):
        path = output / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, separators=(',', ':')) + '\n')
        return {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size}
    files = {}
    for tid, hypotheses in sorted(by_target.items()):
        files[f'targets/{tid}.json'] = write(f'targets/{tid}.json', {'target': targets[tid], 'hypotheses': hypotheses})
    for shard, bundles in sorted(shards.items()):
        files[f'reactions/{shard}.json'] = write(f'reactions/{shard}.json', bundles)
    index = {'schema': 'cannabis-carbon.phase1-hypothesis-view.v1',
        'source_report': 'phase1-target-hypotheses.json', 'source_sha256': hashlib.sha256(raw).hexdigest(),
        'summary': report['summary'], 'claim_boundary': report['claim_boundary'],
        'targets': [{k: t[k] for k in ('cannabisdb_id', 'label', 'status', 'carbon_count', 'structure_status', 'next_step')} |
                    {'hypothesis_count': len(t['hypothesis_ids'])} for t in report['targets']],
        'files': files}
    write('index.json', index)
    print(json.dumps({'targets': len(index['targets']), 'target_lists': len(by_target),
        'reaction_shards': len(shards), 'largest_shard_bytes': max(v['bytes'] for k, v in files.items() if k.startswith('reactions/')),
        'total_bundle_bytes': sum(v['bytes'] for v in files.values())}))


if __name__ == '__main__':
    generate(Path('data/reports/phase1-target-hypotheses.json'), Path('docs/data/hypothesis-view'))
