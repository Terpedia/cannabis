"""Whole-inventory precursor and missing-triglyceride net sensitivity."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build
from .phase1_marts_completions import balanced


def merge_hypotheses(*reports):
    compounds, reactions = {}, {}
    for report in reports:
        for c in report['compounds']:
            if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
                raise ValueError('Conflicting exact compound identity')
            compounds.setdefault(c['id'], c)
        for r in report['reactions']:
            if r['hypothesis_type'] not in ('sn1-acylation', 'sn2-acylation', 'sn3-acylation'):
                raise ValueError('Unexpected precursor hypothesis type')
            if r['id'] in reactions and reactions[r['id']] != r:
                raise ValueError('Conflicting reaction record')
            reactions.setdefault(r['id'], r)
    for r in reactions.values():
        if not balanced([r['left'], r['right']], compounds):
            raise ValueError('Merged precursor equation is not balanced')
    return {'compounds': list(compounds.values()), 'reactions': list(reactions.values())}


def run():
    names = ('full-balanced-network', 'marts-completions', 'catalog-net-gaps',
        'reaction-completion-net', 'lipid-acylation-net', 'triglyceride-net',
        'triglyceride-symmetry-net', 'phosphatidate-hydrolysis-net',
        'glycerolipid-precursors', 'triglyceride-inventory-supplement')
    paths = [Path('data/reports/phase1-' + n + '.json') for n in names]
    inputs = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for doc in inputs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Source snapshot changed')
    network, completions, catalog, completion_net, lipid, tg, symmetry, parent, precursors, supplement = inputs
    original = {**catalog, 'certificates': catalog['certificates'] + [c
                for layer in (completion_net, lipid, tg, symmetry) for c in layer['new_certificates']]}
    report = build(network, completions, original, parent,
                   merge_hypotheses(precursors, supplement), prior_layers=(lipid, tg, symmetry, parent))
    report.update({'schema': 'cannabis-carbon.phase1-glycerolipid-precursors-net.v1',
        'source_sha256': hashes,
        'baseline_certificate_reports': ['phase1-' + n + '.json' for n in names[2:8]],
        'lipid_evidence_reports': ['phase1-glycerolipid-precursors.json', 'phase1-triglyceride-inventory-supplement.json'],
        'claim_boundary': 'Chemistry-only hypotheses, not confirmed Cannabis pathways. '
            'New sn-1, sn-2 and sn-3 acylations are source-forward only; all inherited '
            'direction exclusions and exchange boundaries remain. Existing catalog '
            'and completion equations retain permissive directions. CO2 is the sole '
            'external carbon input. Exact net certificates permit regenerated pre-existing '
            'pools, not established pool startup, energetics, physiological conditions '
            'or compartment feasibility. No enzyme assignment or atom-tracing claim.'})
    report.pop('lipid_evidence_report', None)
    Path('data/reports/phase1-glycerolipid-precursors-net.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
