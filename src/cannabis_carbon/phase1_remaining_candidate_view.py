"""Static all-target search extension, preserving previous scenario bundles."""
import hashlib
import json
from pathlib import Path
from .phase1_thiolase_candidate_view import build as base_view


def build(report, sources):
    bundle = base_view(report, sources, new_candidate_field='new_remaining_gap_candidate')
    bundle['view_scenario'] = 'remaining-candidates'
    bundle['view_boundary'] = '1,609 candidate-linked equations; 161 historical records have conditional net certificates. CO2 is the only net carbon input. Other directions and regenerated pools remain hypotheses, not physiological flux or startup.'
    bundle['restricted_boundary'] = 'Eight reverse exclusions: 153 of 6,220 historical records have conditional net certificates. 5,897 lack producing candidate equations; 135 are solver-infeasible; 35 are exchange species. Candidate proteins are not confirmed enzymes; other directions remain hypothetical.'
    return bundle


def run():
    names = ('phase1-remaining-candidate-net', 'phase1-target-hypotheses', 'phase1-screened-enzyme-overlay',
        'phase1-route-enzyme-overlay', 'phase1-expanded-candidate-net', 'phase1-purine-candidate-net',
        'phase1-replacement-candidate-net', 'phase1-thiolase-candidate-net')
    paths = [Path('data/reports', n + '.json') for n in names]
    reports = [json.loads(p.read_text()) for p in paths]
    for report in reports:
        for p, sha in report.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('View source changed')
    bundle = build(reports[0], reports[1:])
    payload = json.dumps(bundle, separators=(',', ':')) + '\n'
    folder = Path('docs/data/remaining-net-view'); folder.mkdir(parents=True, exist_ok=True)
    (folder / 'bundle.json').write_text(payload)
    manifest = {'schema': 'cannabis-carbon.phase1-remaining-candidate-view.v1', 'file': 'bundle.json',
        'sha256': hashlib.sha256(payload.encode()).hexdigest(), 'bytes': len(payload.encode()),
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'scenario_summaries': {s['id']: s['summary'] for s in reports[0]['scenarios']}}
    (folder / 'index.json').write_text(json.dumps(manifest, separators=(',', ':')) + '\n')
    print(json.dumps({'bytes': manifest['bytes'], 'sha256': manifest['sha256']}))


if __name__ == '__main__':
    run()
