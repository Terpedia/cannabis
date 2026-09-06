"""Full-inventory CO2 net sensitivity for exact PG/PGP reaction proposals."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build
from .phase1_glycerolipid_precursors_net import merge_hypotheses
from .phase1_glycerophospholipid_synthesis import BOUNDARY

TYPES = ('glycerophospholipid-speciation', 'pgp-synthesis', 'pgp-hydrolysis',
         'cdp-dag-synthesis', 'sn1-acylation', 'sn2-acylation')


def run():
    names = ('full-balanced-network', 'marts-completions', 'catalog-net-gaps',
        'reaction-completion-net', 'lipid-acylation-net', 'triglyceride-net',
        'triglyceride-symmetry-net', 'phosphatidate-hydrolysis-net',
        'glycerolipid-precursors-net', 'cardiolipin-net', 'source-mapped-protonation-net',
        'amino-phospholipid-net', 'glycerophospholipid-synthesis')
    paths = [Path('data/reports/phase1-' + n + '.json') for n in names]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p,sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed source snapshot')
    network, completions, catalog = docs[:3]
    original = {**catalog, 'certificates':catalog['certificates'] +
                [c for layer in docs[3:11] for c in layer['new_certificates']]}
    report = build(network, completions, original, docs[11],
        merge_hypotheses(docs[12],allowed_types=TYPES), prior_layers=docs[4:12])
    report.update({'schema':'cannabis-carbon.phase1-glycerophospholipid-net.v1',
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'baseline_certificate_reports':['phase1-' + n + '.json' for n in names[2:12]],
        'lipid_evidence_reports':['phase1-glycerophospholipid-synthesis.json'],
        'claim_boundary':BOUNDARY + ' Full-inventory sensitivity certificates, when present, '
            'establish only exact net conversion under the stated model. All inherited '
            'direction exclusions and external exchanges are retained; CO2 is the sole '
            'external carbon input. Regenerated pre-existing pools are permitted, not '
            'demonstrated zero-pool startup, energetics, compartment feasibility or flux.'})
    report.pop('lipid_evidence_report',None)
    Path('data/reports/phase1-glycerophospholipid-net.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)


if __name__ == '__main__':
    run()
