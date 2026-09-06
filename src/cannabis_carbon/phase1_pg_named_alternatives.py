"""Separate name-derived structural hypotheses; never edits retained PG identities."""
import hashlib
import json
import re
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from .phase1_catalog import stable_id
from .phase1_lipid_acylation import canonical
from .phase1_acyl_input_identity_review import double_stereo_key

BOUNDARY = ('Name-derived structure alternative, not a correction of the historical '
    'target or experimentally resolved identity. Only E/Z geometry explicitly stated '
    'in the CannabisDB XML name is assigned after exact chain-position checks. All '
    'other structure assertions are retained. Original SMILES/InChI omit this geometry '
    'and remain separately preserved. No reaction connects these identity assertions; '
    'no Cannabis detection, pathway or completeness gain is claimed.')


def chains(mol):
    result = {}
    for carbonyl, oxygen, ester, glycerol in mol.GetSubstructMatches(Chem.MolFromSmarts('[C](=[O])-[O]-[C]')):
        position = {2:1,1:2}.get(mol.GetAtomWithIdx(glycerol).GetTotalNumHs())
        if position is None or position in result:
            raise ValueError('Require two distinct terminal/central glycerol ester sites')
        chain = [carbonyl]
        while True:
            options = [a.GetIdx() for a in mol.GetAtomWithIdx(chain[-1]).GetNeighbors()
                       if a.GetAtomicNum()==6 and a.GetIdx() not in chain]
            if not options:
                break
            if len(options)!=1 or mol.GetAtomWithIdx(options[0]).IsInRing():
                raise ValueError('Require unbranched acyclic fatty chain')
            chain.append(options[0])
        result[position] = chain
    if set(result)!={1,2}:
        raise ValueError('Require two glycerol acyl chains')
    return result


def parse_name(name):
    if not name.startswith('PG(') or not name.endswith(')'):
        raise ValueError('Unsupported lipid name')
    parts = name[3:-1].split('/')
    if len(parts)!=2:
        raise ValueError('Require two named chains')
    parsed = []
    for part in parts:
        m = re.fullmatch(r'(\d+):(\d+)(?:\((\d+[EZ](?:,\d+[EZ])*)\))?',part)
        if not m:
            raise ValueError('Unsupported chain notation')
        positions = {} if not m[3] else {int(x[:-1]):x[-1] for x in m[3].split(',')}
        if len(positions)!=int(m[2]):
            raise ValueError('Every double bond requires distinct named geometry')
        parsed.append((int(m[1]),positions))
    return parsed


def derive(mol,name):
    edited = Chem.Mol(mol)
    assignments = []
    for sn,(length,geometry) in enumerate(parse_name(name),1):
        chain = chains(edited)[sn]
        actual = {i+1 for i in range(len(chain)-1)
                  if edited.GetBondBetweenAtoms(chain[i],chain[i+1]).GetBondType()==Chem.BondType.DOUBLE}
        if len(chain)!=length or actual!=set(geometry):
            raise ValueError('Name and encoded chain length/double-bond positions disagree')
        for pos,label in geometry.items():
            if not 2<=pos<len(chain)-1:
                raise ValueError('Require internal alkene')
            bond = edited.GetBondBetweenAtoms(chain[pos-1],chain[pos])
            expected = Chem.BondStereo.STEREOZ if label=='Z' else Chem.BondStereo.STEREOE
            if bond.GetStereo() not in (Chem.BondStereo.STEREONONE,expected):
                raise ValueError('Name conflicts with encoded geometry')
            neighbors = [chain[pos-2],chain[pos+1]]
            if bond.GetBeginAtomIdx()!=chain[pos-1]:
                neighbors.reverse()
            bond.SetStereoAtoms(*neighbors); bond.SetStereo(expected)
            assignments.append({'sn_position':sn,'double_bond_position':pos,'geometry':label})
    Chem.SetDoubleBondNeighborDirections(edited)
    result = Chem.MolFromSmiles(canonical(edited))
    if double_stereo_key(result)!=double_stereo_key(mol):
        raise ValueError('Changed connectivity, charge, isotope or tetrahedral stereo')
    return result,assignments


def run():
    path = Path('data/reports/phase1-pg-source-identity-audit.json')
    audit = json.loads(path.read_bytes())
    for p,sha in audit['source_sha256'].items():
        if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
            raise ValueError('Changed source snapshot')
    rows = []
    for t in audit['targets']:
        source = t['source_fields']
        mol,assignments = derive(Chem.MolFromSmiles(source['smiles']),source['name'])
        smiles = canonical(mol)
        rows.append({'cannabisdb_id':t['cannabisdb_id'],'original_compound_id':t['compound_id'],
            'original_smiles':source['smiles'],'source_name':source['name'],'source_url':t['source_url'],
            'alternative_compound_id':stable_id('structure',smiles),'alternative_smiles':smiles,
            'formula':rdMolDescriptors.CalcMolFormula(mol),'formal_charge':Chem.GetFormalCharge(mol),
            'geometry_assignments':assignments,'identity_status':'name-derived-alternative; original-identity-retained',
            'occurrence_assertion':t['occurrence_assertion'],'claim_boundary':BOUNDARY})
    report = {'schema':'cannabis-carbon.phase1-pg-named-alternatives.v1','targets':rows,'claim_boundary':BOUNDARY,
        'source_sha256':{str(path):hashlib.sha256(path.read_bytes()).hexdigest()},
        'summary':{'source_records':len(rows),'distinct_alternative_structures':len({r['alternative_compound_id'] for r in rows}),
            'changed_historical_identities':0,'new_reactions':0,'new_CO2_route_claims':0}}
    Path('data/reports/phase1-pg-named-alternatives.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)


if __name__ == '__main__':
    run()
