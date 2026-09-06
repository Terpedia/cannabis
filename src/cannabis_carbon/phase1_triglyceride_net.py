"""Full-inventory DGAT sensitivity preserving earlier lipid direction exclusions."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build


def run():
    names = ('phase1-full-balanced-network', 'phase1-marts-completions',
             'phase1-catalog-net-gaps', 'phase1-reaction-completion-net',
             'phase1-lipid-acylation-net', 'phase1-triglyceride-acylation')
    paths = [Path('data/reports', n + '.json') for n in names]
    inputs = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for document in inputs:
        for p, sha in document.get('source_sha256', {}).items():
            if p in hashes and hashes[p] != sha:
                raise ValueError('Changed source snapshot')
    network, completions, catalog, completion_net, lipid_net, hypotheses = inputs
    original = {**catalog, 'certificates': catalog['certificates'] + completion_net['new_certificates']}
    report = build(network, completions, original, lipid_net, hypotheses, prior_layers=(lipid_net,))
    report.update({'schema': 'cannabis-carbon.phase1-triglyceride-net.v1', 'source_sha256': hashes,
        'baseline_certificate_reports': ['phase1-catalog-net-gaps.json', 'phase1-reaction-completion-net.json', 'phase1-lipid-acylation-net.json'],
        'lipid_evidence_report': 'phase1-triglyceride-acylation.json'})
    Path('data/reports/phase1-triglyceride-net.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
