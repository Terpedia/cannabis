"""Find source-backed full equations for incomplete candidate transformations."""
import hashlib
import json
import os
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from .balance import _reaction_smiles_balance


def carbon_participants(smiles):
    """Exact carbon-containing identities; keep stereo, isotopes and charge."""
    if not smiles or smiles.count('>>') != 1 or '*' in smiles:
        return None
    sides = []
    for side in smiles.split('>>'):
        counts = Counter()
        for fragment in side.split('.'):
            molecule = Chem.MolFromSmiles(fragment)
            if molecule is None or molecule.GetNumAtoms() == 0:
                return None
            if any(atom.GetAtomicNum() == 6 for atom in molecule.GetAtoms()):
                for atom in molecule.GetAtoms():
                    atom.SetAtomMapNum(0)
                counts[Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)] += 1
        if not counts:
            return None
        sides.append(dict(sorted(counts.items())))
    return sides


def species_key(sides):
    return tuple(tuple(sorted(side)) for side in sides)


def concrete_participants(smiles):
    sides = []
    for side in smiles.split('>>'):
        counts, details = Counter(), {}
        for fragment in side.split('.'):
            mol = Chem.MolFromSmiles(fragment)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
            identity = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            counts[identity] += 1
            details[identity] = {'smiles': identity, 'formula': rdMolDescriptors.CalcMolFormula(mol),
                                 'formal_charge': Chem.GetFormalCharge(mol),
                                 'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())}
        sides.append([{**details[identity], 'coefficient': coefficient} for identity, coefficient in sorted(counts.items())])
    return sides


def participant_changes(original, reference):
    changes = []
    for side, before, after in zip(('reactants', 'products'), original, reference):
        left, right = {p['smiles']: p for p in before}, {p['smiles']: p for p in after}
        for identity in sorted(left.keys() | right.keys()):
            delta = right.get(identity, {}).get('coefficient', 0) - left.get(identity, {}).get('coefficient', 0)
            if delta:
                changes.append({'side': side, 'coefficient_delta': delta, **{k: v for k, v in (right.get(identity) or left[identity]).items() if k != 'coefficient'}})
    return changes


def match_references(audit, catalog):
    gaps = [r for r in audit['reactions'] if r['status'] != 'balanced']
    keys = {species_key(parts) for r in gaps if (parts := carbon_participants(r['reaction_smarts']))}
    references = {}
    for source in catalog:
        parts = carbon_participants(source['reaction_smarts'])
        if not parts or species_key(parts) not in keys:
            continue
        element, charge = _reaction_smiles_balance(source['reaction_smarts'])
        if not element or not charge or element['status'] != 'balanced' or charge['status'] != 'balanced':
            continue
        references.setdefault(species_key(parts), []).append({**source,
            'carbon_participants': parts, 'participants': concrete_participants(source['reaction_smarts']),
            'element_balance': element, 'charge_balance': charge})
    rows = []
    for gap in gaps:
        original = carbon_participants(gap['reaction_smarts'])
        participants = concrete_participants(gap['reaction_smarts']) if original else None
        matches = references.get(species_key(original), []) if original else []
        candidates = [{**source, 'carbon_stoichiometry_matches': source['carbon_participants'] == original,
                       'participant_changes': participant_changes(participants, source['participants']),
                       'join_method': 'exact-carbon-containing-species-on-each-side; stereo-isotopes-charge-preserved',
                       'required_review': 'Verify omitted carbon-free participants and any changed coefficients against both source records; independently validate enzyme association and physiological direction.'}
                      for source in matches]
        rows.append({'reaction_id': gap['reaction_id'], 'original_reaction_smarts': gap['reaction_smarts'],
                     'original_balance_status': gap['status'], 'original_carbon_participants': original,
                     'original_participants': participants,
                     'balanced_reference_candidates': candidates,
                     'status': 'source-backed-balanced-alternative-found' if candidates else 'no-exact-concrete-reference-found',
                     'claim_boundary': 'Balanced alternatives are separate source equations, not an automatic correction or confirmation of the original MARTS/Rhea transformation. All carbon-free participants and changed coefficients remain explicit in the reference equation.'})
    return rows


def build_report(audit_path, raw, table, query):
    data = raw.read_bytes()
    catalog = json.loads(data)
    print(f'Auditing {len(catalog)} Rhea source equations against balance gaps', flush=True)
    rows = match_references(json.loads(audit_path.read_text()), catalog)
    report = {'schema': 'cannabis-carbon.phase1-balance-reference.v1',
              'generated_at': datetime.now(timezone.utc).isoformat(),
              'source_audit': str(audit_path), 'source_audit_sha256': hashlib.sha256(audit_path.read_bytes()).hexdigest(),
              'catalog': {'table': table, 'query': query, 'snapshot': str(raw), 'sha256': hashlib.sha256(data).hexdigest(), 'row_count': len(catalog)},
              'summary': {'gap_variants': len(rows),
                          'variants_with_balanced_reference_alternatives': sum(bool(r['balanced_reference_candidates']) for r in rows),
                          'variants_with_matching_carbon_stoichiometry': sum(any(c['carbon_stoichiometry_matches'] for c in r['balanced_reference_candidates']) for r in rows),
                          'reference_alternatives': sum(len(r['balanced_reference_candidates']) for r in rows)},
              'rows': rows}
    Path('data/reports/phase1-balance-reference.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(report['summary'])


def run():
    RDLogger.DisableLog('rdApp.warning')
    RDLogger.DisableLog('rdApp.error')
    audit_path = Path('data/reports/terpene-identity-set-candidate-expansion-balance-audit.json')
    table = 'terpedia-489015.terpedia_core.terpene_reaction_smarts_catalog_normalized_current_v2'
    query = f"SELECT rule_id, reaction_smarts, source_url, source_download_url, direction_mode, source_evidence_type FROM `{table}` WHERE STARTS_WITH(rule_id, 'RHEA:')"
    data = subprocess.check_output([os.environ.get('CANNABIS_BQ', 'bq'), '--format=json', 'query',
                                   '--use_legacy_sql=false', '--maximum_bytes_billed=1073741824', '--max_rows=100000', query])
    raw = Path('data/raw/phase1-balance-reference-catalog.json')
    raw.write_bytes(data)
    build_report(audit_path, raw, table, query)


def export_table(report_path, output):
    data = report_path.read_bytes()
    report = json.loads(data)
    records = []
    for row in report['rows']:
        records.append({'reaction_id': row['reaction_id'], 'original_reaction_smarts': row['original_reaction_smarts'],
                        'original_balance_status': row['original_balance_status'], 'review_status': row['status'],
                        'balanced_alternative_count': len(row['balanced_reference_candidates']),
                        'matching_carbon_stoichiometry_count': sum(c['carbon_stoichiometry_matches'] for c in row['balanced_reference_candidates']),
                        'reference_candidates_json': json.dumps(row['balanced_reference_candidates']),
                        'original_participants_json': json.dumps(row['original_participants']),
                        'report_sha256': hashlib.sha256(data).hexdigest(),
                        'catalog_sha256': report['catalog']['sha256'], 'source_audit_sha256': report['source_audit_sha256'],
                        'claim_boundary': row['claim_boundary']})
    output.write_text(''.join(json.dumps(row, separators=(',', ':')) + '\n' for row in records))
    return len(records)


if __name__ == '__main__':
    run()
