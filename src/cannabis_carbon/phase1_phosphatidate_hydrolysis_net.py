"""Full-inventory net sensitivity with source-forward phosphatidate hydrolysis."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build


def run():
    names = ('phase1-full-balanced-network', 'phase1-marts-completions',
        'phase1-catalog-net-gaps', 'phase1-reaction-completion-net',
        'phase1-lipid-acylation-net', 'phase1-triglyceride-net',
        'phase1-triglyceride-symmetry-net', 'phase1-phosphatidate-hydrolysis')
    paths = [Path('data/reports', n + '.json') for n in names]
    inputs = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for doc in inputs:
        for p, sha in doc.get('source_sha256', {}).items():
            if p in hashes and hashes[p] != sha:
                raise ValueError('Source snapshot changed')
    network, completions, catalog, completion_net, lipid_net, tg_net, parent, hypotheses = inputs
    original = {**catalog, 'certificates': catalog['certificates'] + completion_net['new_certificates']
                + lipid_net['new_certificates'] + tg_net['new_certificates']}
    report = build(network, completions, original, parent, hypotheses,
                   prior_layers=(lipid_net, tg_net, parent))
    report.update({'schema': 'cannabis-carbon.phase1-phosphatidate-hydrolysis-net.v1',
        'source_sha256': hashes,
        'baseline_certificate_reports': [n + '.json' for n in names[2:7]],
        'lipid_evidence_report': 'phase1-phosphatidate-hydrolysis.json',
        'claim_boundary': 'Chemistry-only sensitivity, not confirmed Cannabis pathways. '
            'New phosphatidate-hydrolysis equations are source-forward only. All prior '
            'direction exclusions and exchange boundaries are inherited; existing catalog '
            'and completion equations retain permissive directions. Exact net certificates '
            'allow regenerated pre-existing pools, not zero-pool startup, physiological '
            'conditions, energetics or compartment feasibility. No enzyme or atom-tracing claim.'})
    Path('data/reports/phase1-phosphatidate-hydrolysis-net.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
