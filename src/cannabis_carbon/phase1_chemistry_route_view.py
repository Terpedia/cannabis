"""Whole-inventory static Cytoscape certificates for the reaction-first scenario."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_net_view import build as attach_evidence

SOURCES = ('phase1-lipid-acylation-net', 'phase1-reaction-completion-net',
    'phase1-catalog-net-gaps', 'phase1-full-balanced-network', 'phase1-lipid-acylation',
    'phase1-target-hypotheses', 'phase1-screened-enzyme-overlay', 'phase1-route-enzyme-overlay')


def run(triglycerides=False, symmetry=False, hydrolysis=False):
    paths = [Path('data/reports', n + '.json') for n in SOURCES]
    current, parent, baseline, network, lipid, *evidence = [json.loads(p.read_text()) for p in paths]
    previous = current
    layers = [parent, previous]
    extras = []
    if triglycerides or symmetry or hydrolysis:
        extra_paths = [Path('data/reports', n + '.json') for n in ('phase1-triglyceride-net', 'phase1-triglyceride-acylation')]
        current, extra = [json.loads(p.read_text()) for p in extra_paths]
        paths.extend(extra_paths)
        layers.append(current); extras.append(extra)
    if symmetry or hydrolysis:
        extra_paths = [Path('data/reports', n + '.json') for n in ('phase1-triglyceride-symmetry-net', 'phase1-triglyceride-symmetry')]
        current, extra = [json.loads(p.read_text()) for p in extra_paths]
        paths.extend(extra_paths)
        layers.append(current); extras.append(extra)
    if hydrolysis:
        extra_paths = [Path('data/reports', n + '.json') for n in ('phase1-phosphatidate-hydrolysis-net', 'phase1-phosphatidate-hydrolysis')]
        current, extra = [json.loads(p.read_text()) for p in extra_paths]
        paths.extend(extra_paths)
        layers.append(current); extras.append(extra)
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for doc in (*layers, baseline, lipid, *extras):
        for p, sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale chemistry view input')
    certificates = baseline['certificates'] + [c for layer in layers for c in layer['new_certificates']]
    if len({c['compound_id'] for c in certificates}) != len(certificates):
        raise ValueError('Duplicate exact certificate identity')
    used = {s['reaction_id'] for c in certificates for s in c['steps']}
    added_reactions = [r for layer in layers for r in layer['added_reactions']]
    reactions = {r['id']: r for r in network['reactions'] + added_reactions}
    additions = {r['id'] for r in added_reactions}
    source_by_id = {r['rule_id']: r for r in lipid['source_records']}
    for extra in extras:
        source_by_id[extra['source_record']['rule_id']] = extra['source_record']
    selected = []
    for rid in sorted(used):
        r = reactions[rid]
        item = {**r, 'enzyme_evidence_ids': r.get('enzyme_evidence_ids', []),
                'missing_candidate_evidence': not bool(r.get('enzyme_evidence_ids')),
                'is_route_sensitivity': rid in additions}
        if rid in additions:
            item['hypothesis_assumptions'] = [r.get('claim_boundary', parent['claim_boundary'])]
            item['sources'] = r.get('sources', []) + [{'source_urls': [source_by_id[r['source_reaction_id']]['source_url']]
                if r.get('source_reaction_id') in source_by_id else ['https://terpedia.github.io/cannabis/data/reaction-completion-net.json']}]
        selected.append(item)
    cert_by_id = {c['compound_id']: c for c in certificates}
    targets = [{**t, 'certificate_compound_id': t['compound_id'] if t['compound_id'] in cert_by_id else None,
        'startup_status': 'not established by this net certificate',
        'missing_candidate_reaction_ids': sorted({s['reaction_id'] for s in cert_by_id.get(t['compound_id'], {}).get('steps', [])
            if not reactions[s['reaction_id']].get('enzyme_evidence_ids')})} for t in current['targets']]
    required = {p['compound_id'] for r in selected for side in ('left', 'right') for p in r[side]}
    report = {'schema': 'cannabis-carbon.chemistry-route-view.v1', 'view_scenario': 'reaction-first-chemistry',
        'targets': targets, 'certificates': certificates, 'reactions': selected,
        'compounds': [c for c in current['compounds'] if c['id'] in required],
        'summary': {'target_records': len(targets), 'target_status_counts': dict(Counter(t['net_status'] for t in targets))},
        'view_boundary': 'Reaction-first scenario across all 6,220 CannabisDB records. No enzyme gate. Added completion and lipid hypotheses are highlighted as assumptions; every input and coproduct remains in each full equation.',
        'claim_boundary': current['claim_boundary'], 'source_sha256': hashes}
    bundle = attach_evidence(report, evidence)
    folder = Path('docs/data/phosphatidate-hydrolysis-net-view' if hydrolysis else 'docs/data/triglyceride-symmetry-net-view' if symmetry else 'docs/data/triglyceride-net-view' if triglycerides else 'docs/data/chemistry-net-view'); folder.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(bundle, separators=(',', ':')) + '\n'
    (folder / 'bundle.json').write_text(payload)
    manifest = {'schema': report['schema'], 'file': 'bundle.json', 'bytes': len(payload.encode()),
        'sha256': hashlib.sha256(payload.encode()).hexdigest(), 'source_sha256': hashes}
    (folder / 'index.json').write_text(json.dumps(manifest, separators=(',', ':')) + '\n')
    print(json.dumps({'bytes': manifest['bytes'], 'targets': len(targets), 'certificates': len(certificates), 'reactions': len(selected)}))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--triglycerides', action='store_true')
    parser.add_argument('--symmetry', action='store_true')
    parser.add_argument('--hydrolysis', action='store_true')
    args = parser.parse_args()
    run(triglycerides=args.triglycerides, symmetry=args.symmetry, hydrolysis=args.hydrolysis)
