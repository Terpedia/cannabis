"""One shared static bundle for permissive and restricted purine supplements."""
import hashlib
import json
from pathlib import Path
from .phase1_net_view import build as attach_evidence


def build(report, evidence_sources, review):
    scenarios = {s['id']: s for s in report['scenarios']}
    base = {k: v for k, v in report.items() if k not in ('scenarios', 'candidate_reaction_evidence_ids')}
    def scenario(name):
        s = scenarios[name]
        return {**s, 'targets': [{**t, 'startup_status': 'not recomputed; internal pool origin unresolved'} for t in s['targets']]}
    base.update(scenario('permissive-directions'))
    bundle = attach_evidence(base, evidence_sources + [report])
    bundle['restricted_scenario'] = scenario('five-reverse-steps-forbidden')
    bundle['view_scenario'] = 'purine-candidates'
    bundle['view_boundary'] = 'Purine candidate supplement: 109 conditional target certificates. The new 5′-deoxyadenosine route uses unreviewed-reference hypotheses and two previously flagged reverse steps. It disappears under the five-step direction restriction. Neither scenario establishes physiological flux or startup.'
    bundle['restricted_boundary'] = 'Purine candidate supplement with five reverse steps forbidden: 101 conditional target certificates; no gain over the prior restricted scenario. These are analyst-selected restrictions, not established physiological directions. Other directions and regenerated internal pools remain permissive.'
    reviews = {r['reaction_id']: r for r in review['reviews']}
    bundle['reactions'] = [{**r, 'is_new_catalog_candidate': r['new_purine_candidate'],
                            **({'direction_review': reviews[r['id']]} if r['id'] in reviews else {})} for r in bundle['reactions']]
    return bundle


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-purine-candidate-net', 'phase1-target-hypotheses',
        'phase1-screened-enzyme-overlay', 'phase1-route-enzyme-overlay', 'phase1-expanded-candidate-net', 'phase1-candidate-direction-review')]
    reports = [json.loads(p.read_text()) for p in paths]
    for report in (reports[0], reports[-1]):
        for path, sha in report['source_sha256'].items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('View source changed')
    bundle = build(reports[0], reports[1:-1], reports[-1])
    payload = json.dumps(bundle, separators=(',', ':')) + '\n'
    folder = Path('docs/data/purine-net-view'); folder.mkdir(parents=True, exist_ok=True)
    (folder / 'bundle.json').write_text(payload)
    manifest = {'schema': 'cannabis-carbon.phase1-purine-candidate-view.v1', 'file': 'bundle.json',
        'sha256': hashlib.sha256(payload.encode()).hexdigest(), 'bytes': len(payload.encode()),
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'scenario_summaries': {s['id']: s['summary'] for s in reports[0]['scenarios']}}
    (folder / 'index.json').write_text(json.dumps(manifest, separators=(',', ':')) + '\n')
    print(json.dumps({'bytes': manifest['bytes'], 'sha256': manifest['sha256']}))


if __name__ == '__main__':
    run()
