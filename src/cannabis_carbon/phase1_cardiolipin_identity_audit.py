"""Whole-inventory cardiolipin charge and central-stereo review; no reactions added."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from .phase1_lipid_acylation import canonical, exact_scaffold_matches


def build(network, catalog):
    source = next(r for r in catalog if r['rule_id'] == 'RHEA:32932')
    product = Chem.MolFromSmiles(source['reaction_smarts'].split('>>')[1].split('.')[0])
    if Chem.GetFormalCharge(product) != -2:
        raise ValueError('Require source cardiolipin dianion')
    centers = [a.GetIdx() for a in product.GetAtoms() if a.GetAtomicNum() == 6
               and any(n.GetAtomicNum() == 8 and n.GetTotalNumHs() == 1 for n in a.GetNeighbors())]
    if len(centers) != 1:
        raise ValueError('Require unique central glycerol hydroxyl carbon')
    neutral = Chem.RWMol(product)
    for atom in neutral.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
            atom.SetFormalCharge(0); atom.SetNoImplicit(False)
    Chem.SanitizeMol(neutral)
    compounds = {c['id']: c for c in network['compounds']}
    rows = []
    for t in network['targets']:
        mol = Chem.MolFromSmiles(compounds[t['compound_id']]['smiles'])
        matches = exact_scaffold_matches(mol, product)
        state = 'exact-source-dianion'
        if not matches:
            matches = exact_scaffold_matches(mol, neutral)
            state = 'neutral-diacid; explicit-speciation-bridge-required'
        if not matches:
            continue
        stereo = {int(s.centeredOn): str(s.specified) for s in Chem.FindPotentialStereo(mol)}
        statuses = {stereo.get(m[centers[0]], 'achiral') for m in matches}
        if len(statuses) != 1:
            raise ValueError('Ambiguous central-glycerol location')
        rows.append({k: t[k] for k in ('cannabisdb_id', 'compound_id', 'label', 'source_smiles')} | {
            'canonical_smiles': canonical(mol), 'charge_state': state,
            'central_glycerol_stereo': next(iter(statuses)),
            'next_action': 'Construct exact source-compatible PG and CDP-DAG inputs and verify forward reconstruction. Review unresolved central stereo separately; do not infer configuration from names.'})
    return {'schema': 'cannabis-carbon.phase1-cardiolipin-identity-audit.v1',
        'source_record': source, 'targets': rows,
        'claim_boundary': 'Structural review only. Neutral cardiolipin and source dianion are separate identities. '
            'Generic-source product matching does not establish source-compatible precursor stereochemistry, '
            'a balanced exact synthesis reaction, precursor availability or Cannabis activity. No target stereo '
            'has been assigned; no reaction or CO2-route completeness gain is claimed.',
        'summary': {'inventory_target_records_scanned': len(network['targets']),
            'matched_target_records': len(rows),
            'charge_state_counts': dict(Counter(r['charge_state'] for r in rows)),
            'central_glycerol_stereo_counts': dict(Counter(r['central_glycerol_stereo'] for r in rows)),
            'new_reaction_claims': 0, 'new_CO2_route_claims': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-full-balanced-network.json'),
             Path('data/raw/phase1-balance-reference-catalog.json')]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-cardiolipin-identity-audit.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
