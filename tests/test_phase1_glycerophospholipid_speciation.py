import hashlib
import json
from pathlib import Path
from collections import Counter
from rdkit import Chem
from cannabis_carbon.phase1_glycerophospholipid_speciation import build
from cannabis_carbon.phase1_amino_phospholipid_speciation import source_forms
from cannabis_carbon.phase1_lipid_acylation import canonical
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_row_export import encode, decode


def test_full_inventory_audit_reproduces_and_preserves_exact_identities():
    path = Path('data/reports/phase1-glycerophospholipid-speciation.json')
    r = json.loads(path.read_bytes())
    parent = json.loads(Path('data/reports/phase1-amino-phospholipid-net.json').read_bytes())
    catalog = json.loads(Path('data/raw/phase1-balance-reference-catalog.json').read_bytes())
    assert build(parent, catalog) == {k:v for k,v in r.items() if k != 'source_sha256'}
    for p, sha in r['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    assert r['summary']['inventory_records_screened'] == len(parent['targets']) == 6220
    assert r['summary']['new_CO2_route_claims'] == 0
    assert Counter(t['lipid_class'] for t in r['targets']) == {'PG':79, 'PGP':78}
    assert r['summary']['parent_status_counts']['no-net-producing-equation'] == 151
    cs = {c['id']:c for c in r['compounds']}
    original = {t['cannabisdb_id']:t for t in parent['targets']}
    for t in r['targets']:
        assert t['compound_id'] == original[t['cannabisdb_id']]['compound_id']
        assert t['net_status'] == original[t['cannabisdb_id']]['net_status']
        assert t['compound_id'] != t['source_form_compound_id']
        assert balanced([t['left'], t['right']], cs)
        assert t['speciation_type'] == 'explicit-proton-exchange'
        assert t['protons_consumed'] == {'PG':1, 'PGP':3}[t['lipid_class']]
        source = next(s for s in r['source_records'] if s['rule_id'] == t['source_reaction_id'])
        template = next(m for m in map(Chem.MolFromSmiles, source['reaction_smarts'].split('>>')[1].split('.'))
                        if any(a.GetAtomicNum() == 0 for a in m.GetAtoms()))
        mol = Chem.MolFromSmiles(cs[t['compound_id']]['smiles'])
        assert cs[t['source_form_compound_id']]['smiles'] in {canonical(m) for m in source_forms(mol, template)}
        # Removing or inverting encoded headgroup/glycerol stereo must not create a match.
        unassigned = Chem.Mol(mol)
        Chem.RemoveStereochemistry(unassigned)
        assert not source_forms(unassigned, template)
        inverted = Chem.Mol(mol)
        for a in inverted.GetAtoms():
            if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                a.InvertChirality()
        assert not source_forms(inverted, template)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    assert decode(list(reversed(encode(r, sha)))) == r
