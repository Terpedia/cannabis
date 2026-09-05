"""Review-only protonation bridges; never merge identities or alter pathway metrics."""
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize
from .balance import _reaction_smiles_balance
from .phase1_catalog import stable_id

BOUNDARY = ('Structurally inferred proton-transfer hypothesis only. Exact compounds remain distinct; '
    'pH, pKa, compartment, physiological direction and biological relevance are unresolved. '
    'No enzyme is assigned; these records do not change candidate-enzyme or CO2-pathway metrics. '
    'Atom tracing remains deferred.')


def fingerprint(smiles):
    """Lookup key with explicit, independently checked proton-only modifications."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not mol.GetNumAtoms():
        return None, 'invalid-structure'
    if len(Chem.GetMolFrags(mol)) != 1:
        return None, 'multicomponent-excluded'
    if any(a.GetAtomicNum() <= 1 or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None, 'generic-explicit-hydrogen-or-radical-excluded'
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    normalized = rdMolStandardize.Uncharger(True, False, True).uncharge(mol)
    Chem.SanitizeMol(normalized)
    if mol.GetNumAtoms() != normalized.GetNumAtoms():
        return None, 'atom-count-change-excluded'
    atom_signature = lambda a: (a.GetAtomicNum(), a.GetIsotope(), str(a.GetChiralTag()),
                                a.GetIsAromatic(), a.GetNumRadicalElectrons())
    bond_signature = lambda m: [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType()),
        str(b.GetStereo()), tuple(b.GetStereoAtoms()), str(b.GetBondDir())) for b in m.GetBonds()]
    if bond_signature(mol) != bond_signature(normalized):
        return None, 'bond-or-stereo-change-excluded'
    changes = []
    for before, after in zip(mol.GetAtoms(), normalized.GetAtoms()):
        if atom_signature(before) != atom_signature(after):
            return None, 'atom-identity-or-stereo-change-excluded'
        dh = after.GetTotalNumHs() - before.GetTotalNumHs()
        dq = after.GetFormalCharge() - before.GetFormalCharge()
        if dh != dq or ((dh or dq) and before.GetAtomicNum() not in (7, 8, 15, 16)):
            return None, 'non-proton-transfer-change-excluded'
        if dh:
            changes.append({'input_atom_index': before.GetIdx(), 'element': before.GetSymbol(),
                            'hydrogen_delta': dh, 'charge_delta': dq})
    return {'lookup_key': Chem.MolToSmiles(normalized, isomericSmiles=True),
            'canonical_smiles': Chem.MolToSmiles(mol, isomericSmiles=True),
            'formal_charge': Chem.GetFormalCharge(mol),
            'hydrogen_count': sum(a.GetTotalNumHs() for a in mol.GetAtoms()),
            'normalization_changes': changes}, 'eligible'


def bridge(left, right):
    """Return full balanced equation in left-to-right display orientation, not physiology."""
    if left['lookup_key'] != right['lookup_key'] or left['canonical_smiles'] == right['canonical_smiles']:
        return None
    delta = right['formal_charge'] - left['formal_charge']
    if not delta or right['hydrogen_count'] - left['hydrogen_count'] != delta:
        return None  # Deliberately exclude net-zero intramolecular proton relocation.
    sides = [[left['canonical_smiles']], [right['canonical_smiles']]]
    sides[0 if delta > 0 else 1].extend(['[H+]'] * abs(delta))
    equation = '>>'.join('.'.join(side) for side in sides)
    element, charge = _reaction_smiles_balance(equation)
    if not element or not charge or element['status'] != 'balanced' or charge['status'] != 'balanced':
        raise ValueError('Proposed bridge failed independent full-equation balance')
    return {'reaction_smiles': equation, 'protons_consumed': delta,
            'element_balance': element, 'charge_balance': charge}


def build(network):
    compounds = {c['id']: c for c in network['compounds']}
    if len(compounds) != len(network['compounds']):
        raise ValueError('Duplicate compound identifier')
    participants = defaultdict(list)
    for reaction in network['reactions']:
        for side in ('left', 'right'):
            for member in reaction[side]:
                participants[member['compound_id']].append({'reaction_id': reaction['id'],
                    'equation_side': side, 'coefficient': member['coefficient']})
    if not participants.keys() <= compounds.keys():
        raise ValueError('Missing participant identity')
    fingerprints = {cid: fingerprint(c['smiles']) for cid, c in compounds.items()}
    index = defaultdict(list)
    for cid in sorted(participants):
        fp, _ = fingerprints[cid]
        if fp:
            index[fp['lookup_key']].append(cid)
    bridges, targets = {}, []
    for target in network['targets']:
        cid = target['compound_id']
        fp, eligibility = fingerprints[cid]
        bids = []
        if fp:
            for partner in index[fp['lookup_key']]:
                candidate = bridge(fp, fingerprints[partner][0])
                if candidate:
                    bid = stable_id('protonation-bridge', [cid, partner])
                    bids.append(bid)
                    if bid not in bridges:
                        bridges[bid] = {'id': bid, 'target_compound_id': cid,
                            'reaction_participant_compound_id': partner, **candidate,
                            'target_identity_check': fp, 'participant_identity_check': fingerprints[partner][0],
                            'source_reaction_participation': participants[partner],
                            'cannabisdb_ids': [], 'status': 'structurally-inferred-review-required',
                            'enzyme_evidence_ids': [], 'direction_status': 'unresolved; display orientation only',
                            'required_review': ['Check source identity and protonation conventions.',
                                'Measure or source pKa and relevant tissue/compartment pH.',
                                'Determine whether spontaneous equilibration or catalysis applies.',
                                'Only then test a separately labeled all-reactant pathway sensitivity scenario.'],
                            'claim_boundary': BOUNDARY}
                    bridges[bid]['cannabisdb_ids'].append(target['cannabisdb_id'])
        targets.append({k: target[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'source_url')} |
            {'exact_balanced_participation': cid in participants, 'audit_eligibility': eligibility,
             'bridge_ids': bids, 'status': 'protonation-bridge-review-required' if bids else
             'no-protonation-bridge-in-this-audit' if fp else 'excluded-from-protonation-audit'})
    used = {c for b in bridges.values() for c in (b['target_compound_id'], b['reaction_participant_compound_id'])}
    return {'schema': 'cannabis-carbon.phase1-protonation-audit.v1',
        'summary': {'target_records': len(targets), 'bridge_records': len(bridges),
            'targets_with_bridges': sum(bool(t['bridge_ids']) for t in targets),
            'targets_without_exact_participation_with_bridges': sum(bool(t['bridge_ids']) and not t['exact_balanced_participation'] for t in targets),
            'target_status_counts': dict(Counter(t['status'] for t in targets)),
            'structure_eligibility_counts': dict(Counter(status for _, status in fingerprints.values()))},
        'targets': targets, 'bridges': list(bridges.values()),
        'compounds': [compounds[cid] for cid in sorted(used)],
        'method': {'rdkit_version': rdBase.rdkitVersion,
            'uncharger': {'canonicalOrder': True, 'force': False, 'protonationOnly': True},
            'documentation': 'https://www.rdkit.org/docs/cppapi/classRDKit_1_1MolStandardize_1_1Uncharger.html',
            'scope': 'Target-to-balanced-reaction-participant only; net-nonzero proton transfer at N/O/P/S. '
                'No identity merge, tautomer normalization, salt stripping, stereo relaxation, hydride transfer or atom tracing.'},
        'claim_boundary': BOUNDARY}


def run():
    RDLogger.DisableLog('rdApp.warning')
    path = Path('data/reports/phase1-full-balanced-network.json')
    raw = path.read_bytes()
    report = build(json.loads(raw))
    report['source_sha256'] = {str(path): hashlib.sha256(raw).hexdigest()}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    output = Path('data/reports/phase1-protonation-audit.json')
    output.write_text(payload)
    Path('docs/data/phase1-protonation-audit.json').write_text(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    folder = Path('data/derived'); folder.mkdir(exist_ok=True)
    with (folder / 'phase1-protonation-audit.ndjson').open('w') as handle:
        metadata = {k: v for k, v in report.items() if k not in ('targets', 'bridges', 'compounds')}
        for kind, rows in [('metadata', [metadata]), ('target', report['targets']),
                           ('bridge', report['bridges']), ('compound', report['compounds'])]:
            for row in rows:
                handle.write(json.dumps({'record_kind': kind, 'record_id': row.get('id', row.get('cannabisdb_id', 'metadata')),
                    'record_json': json.dumps(row, separators=(',', ':')), 'report_sha256': digest}) + '\n')
    print(json.dumps({'sha256': digest, **report['summary']}))


if __name__ == '__main__':
    run()
