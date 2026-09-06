import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_row_export import encode, decode
from test_phase1_cardiolipin_precursors import hydrolyze, cdp_synthesis, phosphatidyl_transfer
from test_phase1_glycerolipid_precursors import forward_products


def test_every_proposal_replays_and_frontier_is_not_assumed_supplied():
    path = Path('data/reports/phase1-glycerophospholipid-synthesis.json')
    r = json.loads(path.read_bytes())
    for p,sha in r['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    cs = {c['id']:c for c in r['compounds']}
    replay = {'pgp-hydrolysis':hydrolyze, 'pgp-synthesis':phosphatidyl_transfer,
              'cdp-dag-synthesis':cdp_synthesis, 'sn1-acylation':forward_products, 'sn2-acylation':forward_products}
    audit = json.loads(Path('data/reports/phase1-glycerophospholipid-speciation.json').read_bytes())
    bridges = {(json.dumps(t['left'],sort_keys=True),json.dumps(t['right'],sort_keys=True)) for t in audit['targets']}
    for reaction in r['reactions']:
        assert balanced([reaction['left'],reaction['right']],cs)
        assert reaction['enzyme_evidence_ids'] == []
        kind = reaction['hypothesis_type']
        if kind == 'glycerophospholipid-speciation':
            assert (json.dumps(reaction['left'],sort_keys=True),json.dumps(reaction['right'],sort_keys=True)) in bridges
            assert reaction['direction_status'] == 'explicit-reversible-speciation-assumption-not-source-reaction'
        else:
            sides = [[cs[p['compound_id']]['smiles'] for p in reaction[side] for _ in range(p['coefficient'])]
                     for side in ('left','right')]
            products = forward_products(sides[0],kind) if kind.startswith('sn') else replay[kind](sides[0])
            assert products == sorted(sides[1])
            assert reaction['direction_status'] == 'forward-only-generic-source-hypothesis'
    roots = {t['source_form_compound_id'] for t in r['targets']}
    produced = {t['compound_id'] for t in r['precursor_candidates']}
    required = {p for t in r['precursor_candidates'] for p in t['required_precursor_ids']}
    assert roots <= produced
    assert {t['compound_id'] for t in r['frontier']} == (roots | required) - produced
    assert all('not-assumed' in t['status'] for t in r['frontier'])
    assert r['summary']['target_records'] == 157
    assert r['summary']['new_CO2_route_claims'] == 0
    assert decode(list(reversed(encode(r,hashlib.sha256(path.read_bytes()).hexdigest())))) == r
