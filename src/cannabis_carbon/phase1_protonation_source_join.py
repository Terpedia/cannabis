"""Exact endpoint joins to Rhea protonation mappings, without identity promotion."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from rdkit import Chem, RDLogger


def exact_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def build(audit, structures, mappings):
    index = defaultdict(set)
    known = {}
    for chebi, smiles in structures:
        key = exact_key(smiles)
        if chebi in known and known[chebi] != key:
            raise ValueError('Conflicting source structure')
        known[chebi] = key
        if key:
            index[key].add(chebi)
    incoming = defaultdict(list)
    pairs = defaultdict(list)
    for row in mappings:
        if row['ORIGIN'] not in ('curation', 'computation'):
            raise ValueError('Unknown mapping provenance')
        source, major = 'CHEBI:' + row['CHEBI'], 'CHEBI:' + row['CHEBI_PH7_3']
        if not row['CHEBI'].isdigit() or not row['CHEBI_PH7_3'].isdigit():
            raise ValueError('Invalid ChEBI identifier')
        evidence = {'source_chebi_id': source, 'major_chebi_id': major, 'origin': row['ORIGIN']}
        pairs[source, major].append(evidence)
        if source != major:
            incoming[major].append(evidence)
    compounds = {c['id']: c for c in audit['compounds']}
    rows = []
    for bridge in audit['bridges']:
        target = sorted(index.get(exact_key(compounds[bridge['target_compound_id']]['smiles']), set()))
        partner = sorted(index.get(exact_key(compounds[bridge['reaction_participant_compound_id']]['smiles']), set()))
        links = [r for a in target for b in partner for r in pairs[a, b] + pairs[b, a]]
        # Incoming identifiers are retrieval leads only: they do not identify the target.
        leads = [r for b in partner for r in incoming[b] if r['source_chebi_id'] not in known]
        status = ('exact-endpoints-with-source-mapping' if links else
                  'exact-endpoints-without-source-mapping' if target and partner else
                  'missing-exact-target-endpoint' if partner else
                  'missing-exact-participant-endpoint' if target else 'both-exact-endpoints-unresolved')
        rows.append({'id': bridge['id'], 'cannabisdb_ids': bridge['cannabisdb_ids'],
            'target_compound_id': bridge['target_compound_id'],
            'reaction_participant_compound_id': bridge['reaction_participant_compound_id'],
            'target_exact_chebi_ids': target, 'participant_exact_chebi_ids': partner,
            'mapping_evidence': links, 'unresolved_source_structure_leads': leads, 'status': status})
    return {'schema': 'cannabis-carbon.phase1-protonation-source-join.v1', 'rows': rows,
        'summary': {'bridge_records': len(rows), 'status_counts': dict(Counter(r['status'] for r in rows)),
                    'missing_source_structure_ids': len({e['source_chebi_id'] for r in rows
                        for e in r['unresolved_source_structure_leads']})},
        'claim_boundary': 'Exact charged isomeric structures only; no neutralization, tautomer, salt or '
            'stereo relaxation. Missing source endpoint structures remain unresolved. Incoming mapping '
            'identifiers are retrieval leads, not target assignments. Mapping origin is retained; '
            'neither source mapping nor structural agreement establishes Cannabis pH, activity or flux. '
            'No new reaction or net-conversion certificate is asserted.'}


def run():
    RDLogger.DisableLog('rdApp.warning')
    audit_path = Path('data/reports/phase1-expanded-protonation-audit.json')
    receipt_path = Path('data/raw/rhea-protonation-source-retrieval-20260906.json')
    receipt = json.loads(receipt_path.read_text())
    for source in receipt['files']:
        raw = Path(source['path']).read_bytes()
        if len(raw) != source['bytes'] or hashlib.sha256(raw).hexdigest() != source['sha256']:
            raise ValueError('Source snapshot mismatch')
    mapping_path, structure_path = [Path(s['path']) for s in receipt['files']]
    with mapping_path.open() as handle:
        mappings = list(csv.DictReader(handle, delimiter='\t'))
    with structure_path.open() as handle:
        structures = list(csv.reader(handle, delimiter='\t'))
    report = build(json.loads(audit_path.read_text()), structures, mappings)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in (audit_path, receipt_path, mapping_path, structure_path)}
    Path('data/reports/phase1-protonation-source-join.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
