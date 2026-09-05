"""Static candidate-route bundle with preserved synthase and homology evidence."""
import hashlib
import json
from pathlib import Path
from .phase1_net_view import build as attach_evidence


def build(report, sources):
    scenarios = {s['id']: {**s, 'targets': [{**t, 'startup_status': 'not established; regenerated pools allowed'}
        for t in s['targets']]} for s in report['scenarios']}
    links = [{'id': 'synthase-reference-link:' + link['id'],
        'evidence_class': 'reference-linked-candidate; individual sequence evidence retained',
        'screened_proteins': [{'accession': p['candidate_accession'], **p} for p in link['protein_links']],
        'source_report': 'phase1-synthase-reaction-links.json', 'source_link': link}
        for link in report['synthase_reference_links']]
    base = {k: v for k, v in report.items() if k not in ('scenarios', 'candidate_reaction_evidence_ids')}
    base.update(scenarios['permissive-directions'])
    bundle = attach_evidence(base, [*sources, report, {'enzyme_evidence': links}])
    bundle['restricted_scenario'] = scenarios['eight-reverse-steps-forbidden']
    bundle['view_scenario'] = 'thiolase-candidates'
    bundle['view_boundary'] = '1,605 candidate-linked equations; 157 historical target records have conditional net certificates. Only CO2 supplies net carbon. Directions and regenerated pools are hypotheses, not physiological flux or startup.'
    bundle['restricted_boundary'] = 'Eight reverse steps forbidden: 149 of 6,220 historical records have conditional net certificates. 5,897 lack producing candidate equations; 139 are solver-infeasible; 35 are exchange species. Other directions remain hypothetical; not confirmed plant pathways.'
    bundle['reactions'] = [{**r, 'is_new_catalog_candidate': r['new_thiolase_candidate']} for r in bundle['reactions']]
    return bundle


def run():
    names = ('phase1-thiolase-candidate-net', 'phase1-target-hypotheses', 'phase1-screened-enzyme-overlay',
        'phase1-route-enzyme-overlay', 'phase1-expanded-candidate-net', 'phase1-purine-candidate-net',
        'phase1-replacement-candidate-net')
    paths = [Path('data/reports', n + '.json') for n in names]
    reports = [json.loads(p.read_text()) for p in paths]
    for report in reports:
        for p, sha in report.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('View source changed')
    bundle = build(reports[0], reports[1:])
    payload = json.dumps(bundle, separators=(',', ':')) + '\n'
    folder = Path('docs/data/thiolase-net-view'); folder.mkdir(parents=True, exist_ok=True)
    (folder / 'bundle.json').write_text(payload)
    manifest = {'schema': 'cannabis-carbon.phase1-thiolase-candidate-view.v1', 'file': 'bundle.json',
        'sha256': hashlib.sha256(payload.encode()).hexdigest(), 'bytes': len(payload.encode()),
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'scenario_summaries': {s['id']: s['summary'] for s in reports[0]['scenarios']}}
    (folder / 'index.json').write_text(json.dumps(manifest, separators=(',', ':')) + '\n')
    print(json.dumps({'bytes': manifest['bytes'], 'sha256': manifest['sha256']}))


if __name__ == '__main__':
    run()
