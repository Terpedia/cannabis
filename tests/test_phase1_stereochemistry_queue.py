import pytest
import hashlib
import json
from collections import Counter
from pathlib import Path
from cannabis_carbon.phase1_stereochemistry_queue import identity,pair_status
from cannabis_carbon.phase1_stereochemistry_queue import build
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current


def test_stereo_is_not_unspecified_identity_or_protonation():
    assert pair_status('N[C@@H](CO)C(=O)O','N[C@H](CO)C(=O)O')=='specified-stereoisomer-difference'
    assert pair_status('N[C@@H](CO)C(=O)O','NC(CO)C(=O)O')=='unspecified-stereo-identity-review'
    assert pair_status('F/C=C/F','F/C=C\\F')=='specified-stereoisomer-difference'
    assert pair_status('F/C=C/F','FC=CF')=='unspecified-stereo-identity-review'
    assert pair_status('CCO','OCC') is None
    assert pair_status('CC(=O)O','CC(=O)[O-]') is None
    assert pair_status('[13CH3]CO','CCO') is None
    assert pair_status('CC=O','C=CO') is None
    with pytest.raises(ValueError,match='Concrete'):
        identity('*C')


def test_full_model_queue_replay_preserves_equations_identities_and_conflicts():
    root=Path('data/reports'); read=lambda n:json.loads((root/n).read_bytes())
    report=read('phase1-stereochemistry-queue.json'); current=read('phase1-selenium-forward-net.json')
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    layers=[read(n) for n in current['baseline_certificate_reports']]
    reactions,compounds,_=assemble_current(read('phase1-full-balanced-network.json'),read('phase1-marts-completions.json'),current,layers)
    conflicts={r['cannabisdb_id'] for r in read('phase1-no-producer-audit.json')['source_identity_conflicts']}
    rebuilt=build(current,compounds,reactions,conflicts)
    assert rebuilt=={k:v for k,v in report.items() if k!='source_sha256'}
    assert report['summary']['balanced_model_equations']==18177
    assert report['summary']['inventory_records']==6220
    assert report['summary']['coverage_gain_claimed']==report['summary']['new_reactions']==0
    for r in report['reactions']:
        assert r==reactions[r['id']]
    current_targets={t['cannabisdb_id']:t for t in current['targets']}
    for t in report['targets']:
        assert t['net_status']==current_targets[t['cannabisdb_id']]['net_status']=='no-net-producing-equation'
        assert t['source_identity_conflict']==(t['cannabisdb_id'] in conflicts)
        if t['source_identity_conflict']:
            assert t['next_action']=='source-identity-resolution-first'
        for p in t['alternatives']:
            assert p['pair_status']==pair_status(compounds[t['compound_id']]['smiles'],compounds[p['compound_id']]['smiles'])
            assert p['partner_producing_step_ids']
