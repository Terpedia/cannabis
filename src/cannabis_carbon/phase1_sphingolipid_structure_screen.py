"""Inventory structural leads for generic sphingolipid reactions, not substrate assignments."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, rdBase

# Amide-linked amino triol topology; no stereochemical or chain-length assertion.
SMARTS = '[CX3](=[OX1])[NX3][CX4]([CX4][OX2])[CX4]([OX2])[CX4]([OX2])'


def classify(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {'status': 'unparseable', 'match_count': 0}
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(SMARTS), useChirality=False)
    return {'status': 'backbone-topology-lead' if matches else 'no-backbone-match',
            'match_count': len(matches), 'matched_atom_indices': [list(m) for m in matches],
            'unassigned_stereocenters': [i for i, label in Chem.FindMolChiralCenters(mol, includeUnassigned=True) if label == '?']}


def build():
    source = Path('data/reports/phase1-full-balanced-network.json')
    report = json.loads(source.read_text())
    compounds = {c['id']: c for c in report['compounds']}
    results = {cid: classify(compounds[cid]['smiles']) for cid in {t['compound_id'] for t in report['targets']}}
    rows = [{**{k: t[k] for k in ('cannabisdb_id', 'compound_id', 'label', 'structure_status')},
             'smiles': compounds[t['compound_id']]['smiles'], **results[t['compound_id']],
             'model_eligible': False} for t in report['targets']]
    return {'schema': 'cannabis-sphingolipid-structure-screen-v1', 'rows': rows,
        'source_sha256': {str(source): hashlib.sha256(source.read_bytes()).hexdigest()},
        'rdkit_version': rdBase.rdkitVersion, 'smarts': SMARTS, 'use_chirality': False,
        'summary': {'target_records': len(rows), 'exact_structures_screened': len(results),
                    'target_status_counts': dict(Counter(r['status'] for r in rows)),
                    'exact_structure_status_counts': dict(Counter(r['status'] for r in results.values()))},
        'claim_boundary': 'All historical targets retained with unchanged identities. Amide-linked amino-triol topology only: N-acyl chain, base length, C4 stereochemistry, delta-8 position and E/Z geometry are not inferred. Glycosylated or short-chain matches require separate review. No compound is assigned to a generic Rhea substrate from this screen alone; atom tracing remains deferred.'}


if __name__ == '__main__':
    report = build()
    Path('data/reports/phase1-sphingolipid-structure-screen.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))
