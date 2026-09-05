import copy
import hashlib
import json
from collections import Counter
from pathlib import Path
import pytest
from rdkit import Chem
from cannabis_carbon.phase1_catalog import stable_id
from cannabis_carbon.phase1_marts_audit import build


def source(rid, smiles):
    return {'rule_id': rid, 'reaction_smarts': smiles, 'source_uniprot_id': 'REFERENCE_ONLY'}


def test_deduplication_balance_gate_and_exact_identity():
    targets = [{'cannabisdb_id': str(i), 'label': s, 'compound_id': stable_id('structure', s), 'source_url': 'source'}
        for i, s in enumerate(['CCO', 'CC[O-]', 'C[C@H](O)Cl', 'C[C@@H](O)Cl', 'CO'])]
    network = {'targets': targets, 'reactions': []}
    rows = [source('1', 'C=C.O>>CCO'), source('2', 'CCO>>O.C=C'),
            source('3', 'CCO>>CC[O-]'), source('4', 'C[C@H](O)Cl>>C[C@H](O)Cl'), source('5', '*.CO>>CO')]
    result = build(network, rows)
    assert result['summary']['source_balance_status_counts'] == {'balanced': 3, 'imbalanced': 1, 'not-auditable': 1}
    assert result['summary']['balanced_equations'] == 2
    a, b, stereo, other_stereo, generic = result['targets']
    assert a['marts_balanced_matches'] and a['marts_unbalanced_matches']
    assert not b['marts_balanced_matches'] and b['marts_unbalanced_matches']
    assert stereo['marts_balanced_matches'][0]['has_net_production_in_hypothetical_direction'] is False
    assert not other_stereo['marts_balanced_matches']
    assert generic['status'] == 'no-exact-auditable-MARTS-match'
    assert all(r['enzyme_evidence_ids'] == [] for r in result['reactions'])
    with pytest.raises(ValueError, match='Duplicate source'):
        build(network, rows + [rows[0]])


def test_published_whole_marts_audit_and_all_balanced_equations():
    root = Path(__file__).resolve().parents[1]
    raw = (root / 'data/reports/phase1-marts-audit.json').read_bytes()
    assert raw == (root / 'docs/data/phase1-marts-audit.json').read_bytes()
    report = json.loads(raw)
    for name, digest in report['source_sha256'].items():
        path = root / name
        if path.exists():
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    network = json.loads((root / 'data/reports/phase1-full-balanced-network.json').read_text())
    original = copy.deepcopy(network)
    source_rows = [r['source_record'] for r in report['source_ledger']]
    snapshot = root / 'data/raw/phase1-full-marts-snapshot.json'
    if snapshot.exists():
        assert source_rows == json.loads(snapshot.read_text())['rows']
    rebuilt = build(network, source_rows)
    rebuilt['rdkit_version'] = report['rdkit_version']
    assert rebuilt == {k: v for k, v in report.items() if k not in ('source_sha256', 'catalog_provenance')}
    assert network == original
    assert len(source_rows) == 4639 and len(report['targets']) == 6220
    compounds = {c['id']: c for c in report['compounds']}
    for reaction in report['reactions']:
        totals, charges = [], []
        for side in ('left', 'right'):
            atoms, charge = Counter(), 0
            for member in reaction[side]:
                mol = Chem.AddHs(Chem.MolFromSmiles(compounds[member['compound_id']]['smiles']))
                for atom in mol.GetAtoms():
                    atoms[(atom.GetAtomicNum(), atom.GetIsotope())] += member['coefficient']
                charge += Chem.GetFormalCharge(mol) * member['coefficient']
            totals.append(atoms); charges.append(charge)
        assert totals[0] == totals[1] and charges[0] == charges[1]
        assert reaction['enzyme_evidence_ids'] == []
