"""Whole-CannabisDB reaction participation audit without identity relaxation."""
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from rdkit import Chem, RDLogger, rdBase
from .balance import _reaction_smiles_balance


def encoded_structure(smiles):
    if not smiles:
        return None, 'missing-structure'
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or molecule.GetNumAtoms() == 0:
        return None, 'invalid-structure'
    if any(a.GetAtomicNum() == 0 for a in molecule.GetAtoms()):
        return None, 'generic-structure'
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    unspecified = any(str(s.specified) != 'Specified' for s in Chem.FindPotentialStereo(molecule))
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True), (
        'stereo-unspecified-or-unknown' if unspecified else 'no-unassigned-stereo-detected')


def audit_targets(targets, network, catalog, additional_sources=(), matching_ledger_only=False):
    if len({t['id'] for t in targets}) != len(targets):
        raise ValueError('Duplicate CannabisDB source identifier')
    index, reaction_ledger = defaultdict(list), []
    target_structures = [encoded_structure(t.get('smiles')) for t in targets]
    target_keys = {identity for identity, _ in target_structures if identity}
    balance_counts = Counter()
    xref_index = defaultdict(set)
    sources = []
    for reaction in network['reactions']:
        sources.append({'id': 'core:' + reaction['id'], 'source_reaction_id': reaction['id'],
                        'reaction_smiles': reaction.get('reaction_smiles'), 'source_url': reaction.get('source_url'),
                        'source_layer': 'core-networkdb', 'direction_evidence': reaction.get('direction', {}),
                        'orientation_boundary': 'Reaction SMILES orientation; it may differ from separately curated participant direction.'})
        for side in ('reactants', 'products'):
            for member in reaction[side]:
                xref_index[member['compound_id'].lower()].add('core:' + reaction['id'])
    for reaction in catalog['reactions']:
        sources.append({'id': 'expansion:' + reaction['id'], 'source_reaction_id': reaction['source_reaction_id'],
                        'reaction_smiles': reaction['reaction_smarts'], 'source_urls': reaction['source_urls'],
                        'source_layer': 'phase1-expansion-catalog', 'direction_evidence': {},
                        'orientation_boundary': reaction['equation_orientation']})
    sources.extend(additional_sources)
    if len({r['id'] for r in sources}) != len(sources):
        raise ValueError('Duplicate scoped reaction record')
    for source in sources:
        smiles = source['reaction_smiles']
        element, charge = _reaction_smiles_balance(smiles) if smiles and smiles.count('>>') == 1 else (None, None)
        status = 'balanced' if element and charge and element['status'] == charge['status'] == 'balanced' else 'imbalanced' if element and charge else 'not-auditable'
        balance_counts[status] += 1
        participants = defaultdict(list)
        if smiles and smiles.count('>>') == 1:
            for side_name, side in zip(('left', 'right'), smiles.split('>>')):
                counts = Counter()
                for component in side.split('.'):
                    identity, _ = encoded_structure(component)
                    if identity in target_keys:
                        counts[identity] += 1
                for identity, coefficient in counts.items():
                    participants[identity].append({'equation_side': side_name, 'coefficient': coefficient})
        if participants or not matching_ledger_only:
            reaction_ledger.append({**source, 'computed_balance_status': status, 'element_balance': element, 'charge_balance': charge,
                                    'claim_boundary': 'Independently audited source equation; a structure match or balanced equation does not establish physiological direction, enzyme activity, all-input availability or a CO2 pathway.'})
        for identity, roles in participants.items():
            index[identity].append({'reaction_record_id': source['id'], 'source_reaction_id': source['source_reaction_id'],
                                    'computed_balance_status': status, 'roles': roles,
                                    'direction_status': 'not-established-by-structure-match'})
    rows = []
    for target, (identity, structure_status) in zip(targets, target_structures):
        matches = index.get(identity, []) if identity else []
        balanced = [r for r in matches if r['computed_balance_status'] == 'balanced']
        chebi = (target.get('external_ids') or {}).get('chebi')
        chebi_id = 'chebi:' + str(chebi).lower().removeprefix('chebi:') if chebi else None
        xrefs = sorted(xref_index.get(chebi_id, []))
        status = 'balanced-reaction-participant' if balanced else 'reaction-participant-balance-unresolved' if matches else 'no-exact-encoded-reaction-match' if identity else 'structure-unresolved'
        rows.append({'cannabisdb_id': target['id'], 'label': target.get('label'), 'source_url': target.get('source_url'),
                     'source_smiles': target.get('smiles'), 'canonical_isomeric_smiles': identity,
                     'structure_status': structure_status, 'coverage_status': status,
                     'reaction_matches': matches, 'balanced_reaction_record_count': len(balanced),
                     'balanced_right_side_record_count': sum(any(role['equation_side'] == 'right' for role in r['roles']) for r in balanced),
                     'source_chebi_xref': chebi_id, 'xref_reaction_records': xrefs,
                     'xref_boundary': 'CannabisDB-supplied ChEBI cross-reference; not an independently validated structural identity.',
                     'claim_boundary': 'Exact encoded structure only: charge, isotope and stereochemical encoding are retained; unspecified stereo is flagged. Participation and equation side do not prove biosynthesis or CO2 reachability.'})
    return {'schema': 'cannabis-carbon.phase1-target-coverage.v1',
            'summary': {'cannabisdb_records': len(rows), 'coverage_status_counts': dict(Counter(r['coverage_status'] for r in rows)),
                        'structure_status_counts': dict(Counter(r['structure_status'] for r in rows)),
                        'unique_encoded_target_structures': len({r['canonical_isomeric_smiles'] for r in rows if r['canonical_isomeric_smiles']}),
                        'source_reaction_records': len(sources),
                        'reaction_balance_status_counts': dict(balance_counts)},
            'reaction_ledger_scope': 'Only source records with exact target structure participation' if matching_ledger_only else 'All audited source records',
            'targets': rows, 'reaction_ledger': reaction_ledger,
            'metric_scope': 'Every CannabisDB source record is retained. Reaction counts are scoped source records, not deduplicated biochemical reactions. No protonation, tautomer, salt, connectivity-only, or name-based identity merging.',
            'claim_boundary': 'This is reaction-participation coverage, not a complete metabolic map or a claim that carbon comes from CO2. Generic/unbalanced equations and xref-only leads are kept separate; atom tracing remains deferred.'}


def run():
    RDLogger.DisableLog('rdApp.warning')
    RDLogger.DisableLog('rdApp.error')
    paths = [Path('docs/data/compounds.json'), Path('docs/data/networkdb.json'), Path('data/reports/phase1-reaction-catalog.json')]
    targets, network, catalog = [json.loads(path.read_text()) for path in paths]
    result = audit_targets(targets['compounds'], network, catalog)
    result['source_sha256'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    result['rdkit_version'] = rdBase.rdkitVersion
    Path('data/reports/phase1-target-coverage.json').write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(result['summary'])


def export_table(report_path, output):
    data = report_path.read_bytes()
    report = json.loads(data)
    digest = hashlib.sha256(data).hexdigest()
    records = []
    for kind, collection in [('target', 'targets'), ('reaction', 'reaction_ledger')]:
        for row in report[collection]:
            records.append({'record_kind': kind, 'record_id': row.get('cannabisdb_id') or row['id'],
                            'coverage_status': row.get('coverage_status'),
                            'balance_status': row.get('computed_balance_status'),
                            'record_json': json.dumps(row, separators=(',', ':')), 'report_sha256': digest})
    output.write_text(''.join(json.dumps(row, separators=(',', ':')) + '\n' for row in records))
    return len(records)


if __name__ == '__main__':
    run()
