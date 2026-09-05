"""Static adapter for the expanded candidate-only scenario."""
import hashlib
import json
from pathlib import Path
from .phase1_net_view import build


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-expanded-candidate-net',
        'phase1-target-hypotheses', 'phase1-screened-enzyme-overlay', 'phase1-route-enzyme-overlay')]
    reports = [json.loads(p.read_text()) for p in paths]
    for path, sha in reports[0]['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
            raise ValueError('Expanded candidate view source mismatch')
    bundle = build(reports[0], reports[1:] + [reports[0]])
    bundle['view_scenario'] = 'expanded-candidates'
    bundle['view_boundary'] = 'Re-solved candidate-only network: 108 target certificates, seven more than baseline. This is not the frozen full-catalog annotation count.'
    bundle['reactions'] = [{**r, 'is_new_catalog_candidate': r['new_catalog_candidate']} for r in bundle['reactions']]
    # The complete analysis remains in the source report; display only evidence
    # needed by visible equations, not another copy of the whole candidate index.
    for key in ('candidate_reaction_evidence_ids', 'startup_witnesses'):
        bundle.pop(key)
    folder = Path('docs/data/expanded-net-view'); folder.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(bundle, separators=(',', ':')) + '\n'
    (folder / 'bundle.json').write_text(payload)
    manifest = {'schema': 'cannabis-carbon.phase1-expanded-candidate-view.v1', 'file': 'bundle.json',
        'sha256': hashlib.sha256(payload.encode()).hexdigest(), 'bytes': len(payload.encode()),
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': reports[0]['summary']}
    (folder / 'index.json').write_text(json.dumps(manifest, separators=(',', ':')) + '\n')
    print(json.dumps({'bytes': manifest['bytes'], 'sha256': manifest['sha256']}))


if __name__ == '__main__':
    run()
