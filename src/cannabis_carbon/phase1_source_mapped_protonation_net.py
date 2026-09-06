"""Full-inventory protonation sensitivity with inherited chemistry boundaries."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build
from .phase1_source_mapped_protonation import BOUNDARY


def run():
    names = ('full-balanced-network', 'marts-completions', 'catalog-net-gaps',
             'reaction-completion-net', 'lipid-acylation-net', 'triglyceride-net',
             'triglyceride-symmetry-net', 'phosphatidate-hydrolysis-net',
             'glycerolipid-precursors-net', 'cardiolipin-net', 'source-mapped-protonation')
    paths = [Path('data/reports/phase1-' + n + '.json') for n in names]
    docs = [json.loads(p.read_text()) for p in paths]
    for doc in docs:
        for source, sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(source).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed source snapshot')
    network, completions, catalog, completion_net, lipid, tg, symmetry, hydrolysis, precursors, parent, hypotheses = docs
    original = {**catalog, 'certificates': catalog['certificates'] + [c
                for layer in (completion_net, lipid, tg, symmetry, hydrolysis, precursors)
                for c in layer['new_certificates']]}
    report = build(network, completions, original, parent, hypotheses,
                   prior_layers=(lipid, tg, symmetry, hydrolysis, precursors, parent))
    report.update({'schema': 'cannabis-carbon.phase1-source-mapped-protonation-net.v1',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'baseline_certificate_reports': ['phase1-' + n + '.json' for n in names[2:10]],
        'lipid_evidence_reports': ['phase1-source-mapped-protonation.json'],
        'claim_boundary': BOUNDARY + ' This is a separate full-inventory sensitivity scenario. '
            'All inherited direction exclusions and external exchanges remain unchanged. CO2 is '
            'the sole external carbon input. Net certificates permit regenerated pre-existing pools, '
            'not established startup, energetics, transport or physiological flux.'})
    report.pop('lipid_evidence_report', None)
    Path('data/reports/phase1-source-mapped-protonation-net.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
