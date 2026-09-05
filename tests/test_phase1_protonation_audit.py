import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_protonation_audit import fingerprint, bridge, build


def test_proton_transfer_equations_keep_full_stoichiometry():
    acid, _ = fingerprint('CC(=O)O')
    anion, _ = fingerprint('CC(=O)[O-]')
    result = bridge(acid, anion)
    assert result['protons_consumed'] == -1
    assert result['reaction_smiles'] == 'CC(=O)O>>CC(=O)[O-].[H+]'
    assert bridge(anion, acid)['protons_consumed'] == 1
    phosphate, _ = fingerprint('OP(=O)(O)O')
    dianion, _ = fingerprint('OP(=O)([O-])[O-]')
    assert bridge(phosphate, dianion)['protons_consumed'] == -2
    assert bridge(acid, acid) is None


def test_identity_guardrails_exclude_other_chemistry():
    for smiles in ['CC(=O)[O-].[Na+]', '[2H]OC', '[CH3]', '*C']:
        assert fingerprint(smiles)[0] is None
    for a, b in [('C[C@H](O)C(=O)O', 'C[C@@H](O)C(=O)[O-]'),
                 ('[13CH3]C(=O)O', 'CC(=O)[O-]'), ('CC=O', 'C=CO'),
                 ('[CH2+]C', 'CC'), ('C[N+](C)(C)C', 'CN(C)C'),
                 ('C/C=C/C(=O)O', 'C/C=C\\C(=O)[O-]'),
                 ('NCC(=O)O', '[NH3+]CC(=O)[O-]')]:
        aa, _ = fingerprint(a); bb, _ = fingerprint(b)
        assert aa is None or bb is None or bridge(aa, bb) is None


def test_published_audit_replays_all_targets_and_independent_isotope_balance():
    root = Path(__file__).resolve().parents[1]
    path = root / 'data/reports/phase1-protonation-audit.json'
    raw = path.read_bytes(); report = json.loads(raw)
    assert raw == (root / 'docs/data/phase1-protonation-audit.json').read_bytes()
    for name, digest in report['source_sha256'].items():
        source = (root / name).read_bytes()
        assert hashlib.sha256(source).hexdigest() == digest
    network = json.loads(source)
    rebuilt = build(network)
    # Algorithm outcomes must replay; installed RDKit patch version may differ in CI.
    rebuilt['method']['rdkit_version'] = report['method']['rdkit_version']
    assert rebuilt == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert len(report['targets']) == 6220
    assert [t['cannabisdb_id'] for t in report['targets']] == [t['cannabisdb_id'] for t in network['targets']]
    for reaction in report['bridges']:
        totals, charges = [], []
        for side in reaction['reaction_smiles'].split('>>'):
            counts, charge = Counter(), 0
            for fragment in side.split('.'):
                mol = Chem.AddHs(Chem.MolFromSmiles(fragment))
                counts.update((a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms())
                charge += Chem.GetFormalCharge(mol)
            totals.append(counts); charges.append(charge)
        assert totals[0] == totals[1] and charges[0] == charges[1]
        assert reaction['enzyme_evidence_ids'] == []
        assert reaction['status'] == 'structurally-inferred-review-required'
