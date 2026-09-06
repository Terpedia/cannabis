import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_ketone_stereo_hypotheses import ketones,build,REFERENCE
from cannabis_carbon.phase1_marts_completions import balanced


def test_only_defined_secondary_carbon_alcohols_lose_one_stereocentre():
    assert set(ketones('C[C@H](O)CC'))=={'CCC(C)=O'}
    assert ketones('C[C@@H](O)CC').keys()==ketones('C[C@H](O)CC').keys()
    assert not ketones('CC(O)CC')  # no encoded stereocentre
    assert not ketones('CCO')  # primary
    assert not ketones('CC(C)(O)CC')  # tertiary
    assert not ketones('CO[C@H](O)C')  # hemiacetal, not secondary carbon alcohol
    assert not ketones('C[C@H]([O-])CC')  # charged oxygen
    assert not (ketones('C[C@H](O)[C@H](F)CC').keys() & ketones('C[C@@H](O)[C@@H](F)CC').keys())


def test_all_proposals_replay_and_retain_full_cofactor_equations():
    root=Path('data/reports'); read=lambda n:json.loads((root/('phase1-'+n+'.json')).read_bytes())
    report=read('ketone-stereo-hypotheses'); current=read('selenium-forward-net')
    queue=read('stereochemistry-queue'); network=read('full-balanced-network')
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    reference=next(r for r in network['reactions'] if r['id']==REFERENCE)
    assert build(queue,current,reference)=={k:v for k,v in report.items() if k!='source_sha256'}
    compounds={c['id']:c for c in report['compounds']}; reactions={r['id']:r for r in report['reactions']}
    original={c['id'] for c in current['compounds']}
    for r in reactions.values():
        assert balanced([r['left'],r['right']],compounds)
        assert len(r['left'])+len(r['right'])==5
        assert not r['enzyme_evidence_ids']
        assert r['reference_reaction_id']==REFERENCE
    for p in report['pairs']:
        assert p['ketone_already_in_model']==(p['ketone_compound_id'] in original)
        oxid=reactions[p['partner_oxidation_hypothesis_id']]; reduct=reactions[p['target_reduction_hypothesis_id']]
        assert p['partner_compound_id'] in {c['compound_id'] for c in oxid['left']}
        assert p['ketone_compound_id'] in {c['compound_id'] for c in oxid['right']}
        assert p['ketone_compound_id'] in {c['compound_id'] for c in reduct['left']}
        assert p['compound_id'] in {c['compound_id'] for c in reduct['right']}
    assert report['summary']['coverage_gain_claimed']==0
