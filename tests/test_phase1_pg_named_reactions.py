import json
import hashlib
from pathlib import Path
from cannabis_carbon.phase1_marts_completions import balanced
from test_phase1_cardiolipin_precursors import hydrolyze, cdp_synthesis, phosphatidyl_transfer
from test_phase1_glycerolipid_precursors import forward_products


def test_named_reaction_proposals_preserve_originals_and_replay_all_products():
    r=json.loads(Path('data/reports/phase1-pg-named-reactions.json').read_bytes())
    for p,sha in r['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    cs={c['id']:c for c in r['compounds']}
    alts={t['cannabisdb_id']:t for t in r['identity_alternatives']}
    assert len(r['targets'])==23
    for t in r['targets']:
        assert t['compound_id']==alts[t['cannabisdb_id']]['alternative_compound_id']
        assert t['compound_id']!=alts[t['cannabisdb_id']]['original_compound_id']
    original_ids={t['original_compound_id'] for t in alts.values()}
    replay={'pgp-hydrolysis':hydrolyze,'cdp-dag-synthesis':cdp_synthesis,'pgp-synthesis':phosphatidyl_transfer}
    for reaction in r['reactions']:
        assert balanced([reaction['left'],reaction['right']],cs)
        assert not original_ids & {p['compound_id'] for side in ('left','right') for p in reaction[side]}
        kind=reaction['hypothesis_type']
        if kind=='glycerophospholipid-speciation':
            assert reaction['direction_status']=='explicit-reversible-speciation-assumption-not-source-reaction'
            continue
        sides=[[cs[p['compound_id']]['smiles'] for p in reaction[side] for _ in range(p['coefficient'])] for side in ('left','right')]
        products=forward_products(sides[0],kind) if kind.startswith('sn') else replay[kind](sides[0])
        assert products==sorted(sides[1])
        assert reaction['direction_status']=='forward-only-generic-source-hypothesis'
    assert len(r['reactions'])==122
    assert r['summary']['new_CO2_route_claims']==0
