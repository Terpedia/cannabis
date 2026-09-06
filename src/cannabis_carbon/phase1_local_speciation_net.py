"""Full historical inventory sensitivity for explicit local proton hypotheses."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build
from .phase1_glycerolipid_precursors_net import merge_hypotheses
from .phase1_local_speciation_hypotheses import BOUNDARY


def run():
    names = ('full-balanced-network', 'marts-completions', 'catalog-net-gaps',
        'reaction-completion-net', 'lipid-acylation-net', 'triglyceride-net',
        'triglyceride-symmetry-net', 'phosphatidate-hydrolysis-net',
        'glycerolipid-precursors-net', 'cardiolipin-net', 'source-mapped-protonation-net',
        'amino-phospholipid-net', 'glycerophospholipid-net', 'local-speciation-hypotheses')
    paths = [Path('data/reports/phase1-'+n+'.json') for n in names]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed source snapshot')
    network, completions, catalog = docs[:3]
    parent = docs[-2]
    original = {**catalog, 'certificates': catalog['certificates'] +
                [c for layer in docs[3:-2] for c in layer['new_certificates']]}
    if len(original['certificates']) + len(parent['new_certificates']) != 2636:
        raise ValueError('Incomplete prior exact certificates')
    report = build(network, completions, original, parent,
        merge_hypotheses(docs[-1], allowed_types=('local-speciation',)), prior_layers=docs[4:-1])
    if set(report['forbidden_step_ids']) != set(parent['forbidden_step_ids']):
        raise ValueError('Local speciation must retain all inherited direction bounds')
    report.update({'schema':'cannabis-carbon.phase1-local-speciation-net.v1',
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'baseline_certificate_reports':['phase1-'+n+'.json' for n in names[2:-1]],
        'lipid_evidence_reports':['phase1-local-speciation-hypotheses.json'],
        'claim_boundary': BOUNDARY.replace('No CO2-route or historical coverage gain is claimed here.',
            'Full-inventory certificates establish only exact net conversion under these assumptions. '
            'CO2 is the sole external carbon input. Regenerated pre-existing pools are permitted; '
            'startup, energetics, compartment feasibility and Cannabis activity remain unestablished.')})
    report.pop('lipid_evidence_report', None)
    Path('data/reports/phase1-local-speciation-net.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
