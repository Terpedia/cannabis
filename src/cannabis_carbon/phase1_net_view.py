"""Build a compact, source-linked static net-conversion graph bundle."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def build(report, evidence_sources):
    evidence = {}
    for source in evidence_sources:
        for e in source['enzyme_evidence']:
            if e['id'] in evidence and evidence[e['id']] != e:
                raise ValueError('Conflicting evidence identifier')
            evidence[e['id']] = e
    needed = {eid for r in report['reactions'] for eid in r['enzyme_evidence_ids']}
    if not needed <= evidence.keys():
        raise ValueError('Missing reaction enzyme evidence')
    labels = defaultdict(list)
    for t in report['targets']:
        if t['label'] not in labels[t['compound_id']]:
            labels[t['compound_id']].append(t['label'])
    return {**report, 'compounds': [{**c, 'labels': labels[c['id']]} for c in report['compounds']],
        'enzyme_evidence': [evidence[eid] for eid in sorted(needed)]}


def run():
    paths = [Path('data/reports', name + '.json') for name in ['phase1-candidate-net-flux',
        'phase1-target-hypotheses', 'phase1-screened-enzyme-overlay', 'phase1-route-enzyme-overlay']]
    reports = [json.loads(p.read_text()) for p in paths]
    bundle = build(reports[0], reports[1:])
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    # The base report already pins its network and candidate-scope inputs;
    # this manifest additionally pins the exact evidence objects displayed.
    folder = Path('docs/data/net-view'); folder.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(bundle, separators=(',', ':')) + '\n'
    (folder / 'bundle.json').write_text(payload)
    manifest = {'schema': 'cannabis-carbon.phase1-net-view.v1', 'file': 'bundle.json',
        'sha256': hashlib.sha256(payload.encode()).hexdigest(), 'bytes': len(payload.encode()),
        'source_sha256': hashes, 'summary': reports[0]['summary']}
    (folder / 'index.json').write_text(json.dumps(manifest, separators=(',', ':')) + '\n')
    print(json.dumps({'bytes': manifest['bytes'], 'evidence_records': len(bundle['enzyme_evidence']), **manifest['summary']}))


if __name__ == '__main__':
    run()
