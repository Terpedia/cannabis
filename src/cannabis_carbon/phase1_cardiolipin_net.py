"""Full-inventory cardiolipin synthesis/precursor sensitivity, after glycerolipids."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build
from .phase1_glycerolipid_precursors_net import merge_hypotheses

TYPES = ('cardiolipin-synthesis', 'cardiolipin-protonation', 'pgp-hydrolysis',
         'pgp-synthesis', 'cdp-dag-synthesis')


def run():
    names = ('full-balanced-network', 'marts-completions', 'catalog-net-gaps',
        'reaction-completion-net', 'lipid-acylation-net', 'triglyceride-net',
        'triglyceride-symmetry-net', 'phosphatidate-hydrolysis-net',
        'glycerolipid-precursors-net', 'cardiolipin-synthesis', 'cardiolipin-precursors')
    paths = [Path('data/reports/phase1-' + n + '.json') for n in names]
    inputs = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for doc in inputs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Source snapshot changed')
    network, completions, catalog, completion_net, lipid, tg, symmetry, hydrolysis, parent, synthesis, precursors = inputs
    original = {**catalog, 'certificates': catalog['certificates'] + [c
                for layer in (completion_net, lipid, tg, symmetry, hydrolysis) for c in layer['new_certificates']]}
    report = build(network, completions, original, parent,
        merge_hypotheses(synthesis, precursors, allowed_types=TYPES),
        prior_layers=(lipid, tg, symmetry, hydrolysis, parent))
    report.update({'schema': 'cannabis-carbon.phase1-cardiolipin-net.v1',
        'source_sha256': hashes,
        'baseline_certificate_reports': ['phase1-' + n + '.json' for n in names[2:9]],
        'lipid_evidence_reports': ['phase1-cardiolipin-synthesis.json', 'phase1-cardiolipin-precursors.json'],
        'claim_boundary': 'Chemistry-only hypotheses, not confirmed Cannabis pathways. New cardiolipin, '
            'PGP and CDP-DAG synthesis and PGP hydrolysis run only in source-forward orientation. '
            'Explicit cardiolipin protonation equilibria may run both ways. Earlier direction restrictions '
            'and the CO2-only carbon exchange boundary are inherited. Existing catalog/completion '
            'equations retain permissive directions. Exact net certificates permit regenerated '
            'pre-existing pools, not established startup, energetics, compartments or enzyme activity.'})
    report.pop('lipid_evidence_report', None)
    Path('data/reports/phase1-cardiolipin-net.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
