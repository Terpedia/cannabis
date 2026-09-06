"""Bind reviewed mixture-substrate evidence to exact model identities, without edges."""
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from .phase1_catalog import stable_id


def run():
    paths = [Path('data/curation/citral-reductase-evidence.json'),
             Path('data/reports/phase1-alcohol-acetates-net.json'),
             Path('data/reports/phase1-acetate-precursor-gaps.json')]
    evidence, current, gaps = [json.loads(p.read_bytes()) for p in paths]
    for doc in (current, gaps):
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale source')
    compounds = {c['id']: c for c in current['compounds']}
    gap = next(t for t in gaps['targets'] if t['cannabisdb_id'] == 'CDB000585')
    candidates = ['CC(C)=CCC[C@H](C)CC=O', 'CC(C)=CCC/C(C)=C/C=O', 'CC(C)=CCC/C(C)=C\\C=O']
    identities = []
    for sm in candidates:
        mol = Chem.MolFromSmiles(sm); canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        cid = stable_id('structure', canonical)
        identities.append({'compound_id': cid, 'smiles': canonical, 'present_in_model': cid in compounds,
            'tetrahedral_CIP_labels': [label for _, label in Chem.FindMolChiralCenters(mol, includeUnassigned=True)],
            'double_bond_stereo': [str(b.GetStereo()) for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]})
    if identities[0]['tetrahedral_CIP_labels'] != ['S'] or identities[0]['compound_id'] not in gap['target_obstruction']['weights']:
        raise ValueError('Wrong exact citronellal gap identity')
    report = {'schema': 'cannabis-carbon.phase1-citral-reductase-review.v1', 'evidence': evidence,
        'identities': identities, 'affected_target_record': 'CDB000585',
        'summary': {'resolved_model_identity_candidates': sum(t['present_in_model'] for t in identities),
                    'new_reactions': 0, 'coverage_change': 0, 'new_Cannabis_enzyme_assignments': 0},
        'claim_boundary': 'Mixture-substrate evidence does not establish either individual E/Z input channel or enantiopure output. '
            'The S target remains exactly distinct from R citronellal; no bound, compound identity, reaction or coverage changes. '
            'Assay product ratios are evidence observations, not universal flux or stoichiometric constraints.',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    Path('data/reports/phase1-citral-reductase-review.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
