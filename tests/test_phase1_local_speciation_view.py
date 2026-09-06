import hashlib
import json
from collections import Counter
from pathlib import Path


def test_every_sharded_certificate_and_target_matches_full_scientific_results():
    root=Path('docs/data/local-speciation-net-view')
    def read(name):
        return json.loads(Path(name).read_bytes())
    def checked(ref):
        path=root/ref['file']
        assert path.resolve().is_relative_to(root.resolve())
        payload=path.read_bytes()
        assert len(payload)==ref['bytes']
        assert hashlib.sha256(payload).hexdigest()==ref['sha256']
        return json.loads(payload)
    index=checked(read(root/'index.json'))
    shared=checked(index['shared_chemistry'])
    current=read('data/reports/phase1-local-speciation-net.json')
    previous=read('docs/data/glycerophospholipid-net-view/bundle.json')
    for p,sha in index['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    expected={c['compound_id']:c for c in previous['certificates']+current['new_certificates']}
    assert len(expected)==len(index['certificate_files'])==2720
    assert set(index['certificate_files'])==set(expected)
    for cid,ref in index['certificate_files'].items():
        assert checked(ref)==expected[cid]
        assert ref['bytes']<50000
    assert (root/'bundle.json').stat().st_size<5000000
    assert index['shared_chemistry']['bytes']<16000000
    assert index['summary']['target_records']==6220
    assert index['model_summary']==current['summary']
    assert index['summary']['target_status_counts']==current['summary']['net_status_counts']
    for k in ('forbidden_step_ids','co2_compound_id','external_exchange_compound_ids','claim_boundary'):
        assert index[k]==current[k]
    reactions={r['id']:r for r in shared['reactions']}
    compounds={c['id']:c for c in shared['compounds']}
    evidence={e['id'] for e in shared['enzyme_evidence']}
    assert len(reactions)==len(shared['reactions'])
    required={s['reaction_id'] for c in expected.values() for s in c['steps']}
    assert set(reactions)==required
    for r in reactions.values():
        assert r['sources']
        assert set(r['enzyme_evidence_ids'])<=evidence
        assert {p['compound_id'] for side in ('left','right') for p in r[side]}<=compounds.keys()
        if r.get('hypothesis_type')=='local-speciation':
            assert r['is_route_sensitivity'] is True
            assert r['claim_boundary'] in r['hypothesis_assumptions']
            assert any(s.get('evidence_type')=='computed-local-proton-change-not-curated-reaction' for s in r['sources'])
    for old in previous['reactions']:
        assert reactions[old['id']]==old
    for original in current['certificate_reactions']:
        r=reactions[original['id']]
        for side in ('left','right'):
            assert r[side]==original[side]
    for target, source in zip(index['targets'], current['targets'],strict=True):
        assert all(target[k]==v for k,v in source.items())
        cid=source['compound_id']
        assert target['certificate_compound_id']==(cid if cid in expected else None)
        gaps={s['reaction_id'] for s in expected.get(cid,{}).get('steps',[]) if not reactions[s['reaction_id']]['enzyme_evidence_ids']}
        assert target['missing_candidate_reaction_count']==len(gaps)
    assert Counter(t['net_status'] for t in index['targets'])['exact-net-conversion-hypothesis']==2723
