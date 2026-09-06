import hashlib
import json
from collections import Counter
from pathlib import Path
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_row_export import encode,decode


def test_every_paired_probe_preserves_identity_boundaries_and_exact_net_balance():
    path=Path('data/reports/phase1-pg-named-net.json')
    assert path.exists(), 'Paired calculation is pending; no result may be published yet'
    r=json.loads(path.read_bytes())
    current=json.loads(Path('data/reports/phase1-glycerophospholipid-net.json').read_bytes())
    proposals=json.loads(Path('data/reports/phase1-pg-named-reactions.json').read_bytes())
    for p,sha in r['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    assert r['historical_inventory_summary_unchanged']==current['summary']
    assert r['summary']['historical_coverage_gain_claimed']==r['summary']['changed_historical_identities']==0
    for k in ('external_exchange_compound_ids','co2_compound_id'):
        assert r[k]==current[k]
    cs={c['id']:c for c in r['compounds']}
    assert {cid for cid in r['external_exchange_compound_ids'] if cs[cid]['carbon_count']}=={r['co2_compound_id']}
    assert cs[r['co2_compound_id']]['smiles']=='O=C=O'
    forbidden=set(r['forbidden_step_ids'])
    assert forbidden==set(current['forbidden_step_ids']) | {
        a['id']+':hypothetical-right-to-left' for a in r['added_reactions']
        if a['hypothesis_type']!='glycerophospholipid-speciation'}
    allreactions={a['id']:a for a in r['certificate_reactions']+r['added_reactions']}
    steps={s['id']:s for s in orientations(list(allreactions.values())) if s['id'] not in forbidden}
    source={a['id']:a for a in proposals['reactions']}
    for reaction in allreactions.values():
        assert balanced([reaction['left'],reaction['right']],cs)
    for a in r['added_reactions']:
        assert a==source[a['id']]
    assert {a['id'] for a in r['added_reactions']} | {j['hypothesis_id'] for j in r['existing_equation_joins']}==set(source)
    alts={t['cannabisdb_id']:t for t in proposals['identity_alternatives']}
    assert len(r['paired_probes'])==len(alts)==23
    original_ids={t['original_compound_id'] for t in alts.values()}
    for a in r['added_reactions']:
        assert not original_ids & {p['compound_id'] for s in ('left','right') for p in a[s]}
    for pair in r['paired_probes']:
        identity=alts[pair['cannabisdb_id']]
        assert pair['original_result']['compound_id']==identity['original_compound_id']
        for k in ('alternative_baseline_result','alternative_extended_result'):
            assert pair[k]['compound_id']==identity['alternative_compound_id']
        assert identity['original_compound_id']!=identity['alternative_compound_id']
        for k in ('original_result','alternative_baseline_result','alternative_extended_result'):
            cert=pair[k]
            if cert['status']=='exact-net-conversion-hypothesis':
                validate_certificate(cert,steps,cs,set(r['external_exchange_compound_ids']),r['co2_compound_id'])
                assert not forbidden & {s['step_id'] for s in cert['steps']}
                if k=='alternative_baseline_result':
                    assert not {s['reaction_id'] for s in cert['steps']} & {a['id'] for a in r['added_reactions']}
    for k in ('original_result','alternative_baseline_result','alternative_extended_result'):
        assert r['summary'][k+'_counts']==dict(Counter(p[k]['status'] for p in r['paired_probes']))
    assert r['summary']['balanced_equations']==current['summary']['balanced_equations']+len(r['added_reactions'])
    assert decode(encode(r,hashlib.sha256(path.read_bytes()).hexdigest())[::-1])==r
