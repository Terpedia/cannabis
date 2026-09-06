"""Full-inventory net sensitivity for exact amino-phospholipid reaction gaps."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build
from .phase1_glycerolipid_precursors_net import merge_hypotheses
from .phase1_amino_phospholipid_synthesis import BOUNDARY

TYPES = ('pe-synthesis', 'ps-synthesis', 'pe-first-methylation', 'pe-second-methylation',
         'amino-phospholipid-speciation', 'cdp-dag-synthesis', 'sn1-acylation', 'sn2-acylation')


def run():
    names = ('full-balanced-network', 'marts-completions', 'catalog-net-gaps',
             'reaction-completion-net', 'lipid-acylation-net', 'triglyceride-net',
             'triglyceride-symmetry-net', 'phosphatidate-hydrolysis-net',
             'glycerolipid-precursors-net', 'cardiolipin-net', 'source-mapped-protonation-net',
             'amino-phospholipid-synthesis', 'amino-phospholipid-precursors')
    paths = [Path('data/reports/phase1-' + n + '.json') for n in names]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p, sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed source snapshot')
    network, completions, catalog = docs[:3]
    parent = docs[10]
    original = {**catalog, 'certificates': catalog['certificates'] +
                [c for layer in docs[3:10] for c in layer['new_certificates']]}
    hypotheses = merge_hypotheses(*docs[11:], allowed_types=TYPES)
    report = build(network, completions, original, parent, hypotheses, prior_layers=docs[4:11])
    report.update({'schema': 'cannabis-carbon.phase1-amino-phospholipid-net.v1',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'baseline_certificate_reports': ['phase1-' + n + '.json' for n in names[2:11]],
        'lipid_evidence_reports': ['phase1-' + n + '.json' for n in names[11:]],
        'claim_boundary': BOUNDARY + ' Full-inventory chemistry sensitivity. New synthesis '
            'and methylation reactions are source-forward only; separate speciation steps '
            'are assumed reversible. All inherited exclusions and external exchanges are '
            'unchanged. CO2 is the sole external carbon input. Certificates permit regenerated '
            'pre-existing pools, not demonstrated zero-pool startup or physiological flux.'})
    report.pop('lipid_evidence_report', None)
    Path('data/reports/phase1-amino-phospholipid-net.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
