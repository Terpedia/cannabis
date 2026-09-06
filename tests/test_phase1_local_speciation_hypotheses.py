import hashlib
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from cannabis_carbon.phase1_catalog import stable_id
from cannabis_carbon.phase1_local_speciation_hypotheses import build, proposal
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_scope import orientations


def compounds_for(smiles):
    records = {}
    ids = []
    for s in smiles + ['[H+]']:
        mol = Chem.MolFromSmiles(s)
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        cid = stable_id('structure', canonical)
        records[cid] = {'id':cid, 'smiles':canonical, 'formal_charge':Chem.GetFormalCharge(mol),
                        'formula':rdMolDescriptors.CalcMolFormula(mol)}
        ids.append(cid)
    return ids, records


def test_local_speciation_keeps_stereo_isotopes_bonds_and_explicit_proton_balance():
    for left, right, delta in [('NCC(=O)O','[NH3+]CC(=O)[O-]',0),
                                ('CN','C[NH3+]',1), ('CC(=O)O','CC(=O)[O-]',-1)]:
        ids, cs = compounds_for([left,right])
        r = proposal(ids[0],ids[1],cs)
        assert r and r['protons_consumed'] == delta
        assert balanced([r['left'],r['right']],cs)
        assert len(r['left'])+len(r['right']) == (3 if delta else 2)
    for left, right in [('N[C@@H](C)C(=O)O','[NH3+][C@H](C)C(=O)[O-]'),
                        ('NCC(=O)O','[NH3+][13CH2]C(=O)[O-]'),
                        ('CC=O','C=CO'), ('CN','CN'), ('C/C=C/CN','C/C=C\\C[NH3+]')]:
        ids, cs = compounds_for([left,right])
        assert proposal(ids[0],ids[1],cs) is None


def test_all_speciation_hypotheses_replay_and_retain_partner_provenance():
    path = Path('data/reports/phase1-local-speciation-hypotheses.json')
    report = json.loads(path.read_bytes())
    gaps = json.loads(Path('data/reports/phase1-current-reaction-gaps.json').read_bytes())
    assert build(gaps) == {k:v for k,v in report.items() if k!='source_sha256'}
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    assert report['summary']['balanced_hypotheses'] == 92
    assert report['summary']['target_records'] == 93
    assert report['summary']['speciation_type_counts'] == {
        'net-zero-intramolecular-proton-relocation':50, 'explicit-proton-exchange':42}
    assert report['summary']['coverage_gain_claimed'] == 0
    cs = {c['id']:c for c in report['compounds']}
    source = {s['id']:s for s in orientations(report['partner_source_reactions'])}
    conflicts = {t['cannabisdb_id'] for t in gaps['targets'] if t['category']=='source-identity-conflict'}
    for r in report['reactions']:
        assert not conflicts & set(r['cannabisdb_ids'])
        assert balanced([r['left'],r['right']],cs)
        assert r['target_compound_id'] != r['reaction_participant_compound_id']
        assert not r['enzyme_evidence_ids']
        assert r['source_evidence_type']=='computed-local-proton-change-not-curated-reaction'
        for sid in r['partner_producing_step_ids']:
            assert sid in source and sid not in gaps['forbidden_step_ids']
            step = source[sid]; cid = r['reaction_participant_compound_id']
            assert sum(p['coefficient'] for p in step['outputs'] if p['compound_id']==cid) > sum(
                p['coefficient'] for p in step['required_inputs'] if p['compound_id']==cid)
