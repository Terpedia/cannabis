"""Review unresolved acyl inputs without converting identity ambiguity into reactions."""
import hashlib
import json
from collections import defaultdict, Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from .phase1_lipid_acylation import canonical

BOUNDARY = ('Identity-review leads only. Comparison removes only double-bond stereo '
    'for indexing, retaining connectivity, isotopes, charge and tetrahedral stereo. '
    'A comparison hit is not an exact identity match, reaction, isomerization, '
    'protonation bridge or permission to assign missing E/Z geometry. Candidate '
    'structures remain distinct; no target identity or pathway coverage is changed.')


def double_stereo_key(mol):
    copy = Chem.Mol(mol)
    for bond in copy.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            bond.SetStereo(Chem.BondStereo.STEREONONE)
        if bond.GetBondDir() in (Chem.BondDir.ENDUPRIGHT,Chem.BondDir.ENDDOWNRIGHT):
            bond.SetBondDir(Chem.BondDir.NONE)
    return canonical(copy)


def alkene_geometry(mol):
    return [{'bond_index':b.GetIdx(),'atom_indices':[b.GetBeginAtomIdx(),b.GetEndAtomIdx()],
             'stereo':str(b.GetStereo())} for b in mol.GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE and
            b.GetBeginAtom().GetAtomicNum() == b.GetEndAtom().GetAtomicNum() == 6]


def build(audit, current):
    cs = {c['id']:c for c in current['compounds']}
    unresolved = [r['compound_id'] for r in audit['input_results']
                  if r['status'] == 'no-net-producing-candidate-equation']
    keys = {cid:double_stereo_key(Chem.MolFromSmiles(cs[cid]['smiles'])) for cid in unresolved}
    index = defaultdict(list)
    for c in current['compounds']:
        mol = Chem.MolFromSmiles(c['smiles'])
        key = double_stereo_key(mol)
        if key in set(keys.values()):
            index[key].append(c['id'])
    input_status = {r['compound_id']:r['status'] for r in audit['input_results']}
    rows = []
    for cid in unresolved:
        geometry = alkene_geometry(Chem.MolFromSmiles(cs[cid]['smiles']))
        candidates = []
        for other in sorted(index[keys[cid]]):
            if other == cid:
                continue
            candidates.append({'compound_id':other,'smiles':cs[other]['smiles'],
                'alkene_geometry':alkene_geometry(Chem.MolFromSmiles(cs[other]['smiles'])),
                'individual_supply_status':input_status.get(other,'not-tested-by-input-audit'),
                'identity_relationship':'double-bond-stereo-ignored-comparison; not-exact-identity'})
        rows.append({'compound_id':cid,'smiles':cs[cid]['smiles'],'alkene_geometry':geometry,
            'unassigned_alkene_bonds':sum(g['stereo']=='STEREONONE' for g in geometry),
            'comparison_key':keys[cid],'candidate_structures':candidates,
            'affected_cannabisdb_ids':[t['cannabisdb_id'] for t in audit['targets'] if cid in t['unresolved_terminal_input_ids']],
            'next_curation_step':'Retrieve source structure and supporting characterization for affected CannabisDB records; resolve E/Z geometry or retain explicitly ambiguous identity. Do not add a reaction from a defined isomer to an unspecified structure.',
            'claim_boundary':BOUNDARY})
    return {'schema':'cannabis-carbon.phase1-acyl-input-identity-review.v1','claim_boundary':BOUNDARY,
        'inputs':rows,'summary':{'inventory_structures_screened':len(cs),'unresolved_input_structures':len(rows),
            'inputs_with_unassigned_alkene_geometry':sum(r['unassigned_alkene_bonds']>0 for r in rows),
            'comparison_candidate_counts':dict(Counter(str(len(r['candidate_structures'])) for r in rows)),
            'affected_target_records':len({t for r in rows for t in r['affected_cannabisdb_ids']}),
            'new_reactions':0,'new_CO2_route_claims':0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-glycerophospholipid-input-audit.json'),
             Path('data/reports/phase1-glycerophospholipid-net.json')]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p,sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed input snapshot')
    report = build(*docs)
    report['source_sha256'] = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-acyl-input-identity-review.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)


if __name__ == '__main__':
    run()
