"""Exact primary/secondary monoacetate candidates; no inferred enzyme assignment."""
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced

REFERENCE = 'balanced-equation:b2ec3259c31c217eaad784236834d1d76101c7fa619d5694b20cd6546ceb772c'
ACETATE = Chem.MolFromSmarts('[CH3][CX3](=[OX1])[OX2][#6]')
ESTER = Chem.MolFromSmarts('[CX3](=[OX1])[OX2][#6]')
BOUNDARY = ('Primary/secondary alcohol monoacetate reaction-class hypothesis. Exact deacetylated '
    'alcohol identity, stereochemistry and acetyl-CoA/CoA charged structures are retained. '
    'Terpedia benzyl-alcohol acetyltransfer chemistry is a template, not evidence of the proposed '
    'substrate scope or Cannabis activity. All new equations are proposed forward-only. '
    'No enzyme assignment or full CO2 pathway follows from this proposal alone.')


def precursor(smiles):
    mol = Chem.MolFromSmiles(smiles)
    matches = mol.GetSubstructMatches(ACETATE)
    if not matches:
        return None, 'not-an-acetate-ester'
    if len(matches) != 1 or len(mol.GetSubstructMatches(ESTER)) != 1:
        return None, 'multiple-esters-separate-site-specific-review'
    methyl, carbonyl, oxygen, bridge, alcohol_carbon = matches[0]
    carbon = mol.GetAtomWithIdx(alcohol_carbon)
    if carbon.GetIsAromatic() or carbon.GetHybridization() != Chem.HybridizationType.SP3:
        return None, 'phenolic-or-unsaturated-oxygen-substrate-separate-review'
    if sum(a.GetAtomicNum() == 6 for a in carbon.GetNeighbors()) > 2:
        return None, 'tertiary-alcohol-substrate-separate-review'
    if Chem.GetFormalCharge(mol) or any(a.GetIsotope() for a in mol.GetAtoms()):
        return None, 'charged-or-isotope-specific-target-separate-review'
    edit = Chem.RWMol(mol)
    for idx in sorted((methyl, carbonyl, oxygen), reverse=True):
        edit.RemoveAtom(idx)
    product = edit.GetMol(); Chem.SanitizeMol(product)
    if len(Chem.GetMolFrags(product)) != 1:
        raise ValueError('Unexpected disconnected alcohol')
    return Chem.MolToSmiles(product, isomericSmiles=True), 'exact-deacetylated-alcohol'


def build(current, network):
    compounds = {c['id']: c for c in current['compounds']}
    ref = next(r for r in network['reactions'] if r['id'] == REFERENCE)
    alcohol = next(p['compound_id'] for p in ref['left'] if compounds[p['compound_id']]['smiles'] == 'OCc1ccccc1')
    ester = next(p['compound_id'] for p in ref['right'] if compounds[p['compound_id']]['smiles'] == 'CC(=O)OCc1ccccc1')
    rows = []; reactions = {}; used = set()
    for t in current['targets']:
        if t['net_status'] != 'no-net-producing-equation':
            continue
        sm = compounds[t['compound_id']]['smiles']
        derived, status = precursor(sm)
        if status == 'not-an-acetate-ester':
            continue
        row = {k: t[k] for k in ('cannabisdb_id', 'label', 'compound_id')}
        row.update({'target_smiles': sm, 'status': status})
        if derived:
            cid = stable_id('structure', derived)
            row.update({'alcohol_smiles': derived, 'alcohol_compound_id': cid, 'alcohol_present_in_model': cid in compounds})
            if cid not in compounds:
                row['status'] = 'exact-alcohol-absent-from-model'
            else:
                mapping = {alcohol: cid, ester: t['compound_id']}
                sides = [[{**p, 'compound_id': mapping.get(p['compound_id'], p['compound_id'])} for p in ref[s]] for s in ('left', 'right')]
                if not balanced(sides, compounds):
                    raise ValueError('Unbalanced acetate proposal')
                rid = stable_id('alcohol-acetate-hypothesis', sides)
                reactions[rid] = {'id': rid, 'left': sides[0], 'right': sides[1],
                    'hypothesis_type': 'alcohol-acetylation', 'reference_reaction_id': REFERENCE,
                    'source_url': ref['sources'][0]['source_urls'][0],
                    'source_evidence_type': 'reaction-class-analogy-not-target-substrate-or-Cannabis-evidence',
                    'direction_status': 'proposed-forward-only', 'enzyme_evidence_ids': [],
                    'balance_status': 'independently-element-isotope-charge-balanced', 'claim_boundary': BOUNDARY}
                row['hypothesis_id'] = rid
                row['status'] = 'balanced-proposal-with-existing-exact-precursors'
                used.update(p['compound_id'] for s in sides for p in s)
        rows.append(row)
    return {'schema': 'cannabis-carbon.phase1-alcohol-acetates.v1', 'targets': rows,
        'reactions': list(reactions.values()), 'reference_reactions': [ref],
        'compounds': [compounds[c] for c in sorted(used)], 'claim_boundary': BOUNDARY,
        'summary': {'acetate_target_records_reviewed': len(rows), 'balanced_proposed_equations': len(reactions),
                    'new_compound_structures': 0, 'coverage_gain_claimed': 0}}


def run():
    paths = [Path('data/reports/phase1-' + n + '.json') for n in ('c17-elongation-net', 'full-balanced-network')]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale source')
    report = build(*docs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-alcohol-acetates.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
