import hashlib
import json
from collections import Counter
from pathlib import Path
from unittest.mock import patch
import pytest
from cannabis_carbon.phase1_glycerophospholipid_net import TYPES
from cannabis_carbon.phase1_glycerolipid_precursors_net import merge_hypotheses
from cannabis_carbon.phase1_lipid_acylation_net import build, equation_key
from cannabis_carbon.phase1_marts_completions import balanced
from cannabis_carbon.phase1_scope import orientations
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_row_export import encode, decode


def test_new_synthesis_forward_only_and_speciation_explicitly_reversible():
    evidence = json.loads(Path('data/reports/phase1-glycerophospholipid-synthesis.json').read_bytes())
    merged = merge_hypotheses(evidence,allowed_types=TYPES)
    compounds = {c['id']:c for c in merged['compounds']}
    compounds['co2'] = {'id':'co2','smiles':'O=C=O','carbon_count':1}
    parent = {'targets':[], 'external_exchange_compound_ids':['co2'], 'co2_compound_id':'co2'}
    prior = {'compounds':[], 'added_reactions':[], 'forbidden_step_ids':['inherited:reverse']}
    prefix = 'cannabis_carbon.phase1_lipid_acylation_net.'
    with patch(prefix+'assemble',return_value=({},compounds,[])), \
         patch(prefix+'NetModel',side_effect=RuntimeError('captured')) as model:
        with pytest.raises(RuntimeError,match='captured'):
            build({'targets':[]},{},parent,parent,merged,prior_layers=(prior,))
    unique = {}
    for r in merged['reactions']:
        unique.setdefault(equation_key(r),r)
    assert set(model.call_args.kwargs['forbidden_step_ids']) == {'inherited:reverse'} | {
        r['id']+':hypothetical-right-to-left' for r in unique.values()
        if r['hypothesis_type'] != 'glycerophospholipid-speciation'}


def test_full_inventory_and_every_new_exact_certificate():
    path = Path('data/reports/phase1-glycerophospholipid-net.json')
    assert path.exists(), 'Full calculation pending; no coverage claim may be published'
    report = json.loads(path.read_bytes())
    parent = json.loads(Path('data/reports/phase1-amino-phospholipid-net.json').read_bytes())
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    assert [(t['cannabisdb_id'],t['compound_id']) for t in report['targets']] == [
        (t['cannabisdb_id'],t['compound_id']) for t in parent['targets']]
    assert len(report['targets']) == len({t['cannabisdb_id'] for t in report['targets']}) == 6220
    for k in ('co2_compound_id','external_exchange_compound_ids'):
        assert report[k] == parent[k]
    cs = {c['id']:c for c in report['compounds']}
    assert {c for c in report['external_exchange_compound_ids'] if cs[c]['carbon_count']} == {report['co2_compound_id']}
    assert cs[report['co2_compound_id']]['smiles'] == 'O=C=O'
    reactions = {r['id']:r for r in report['added_reactions']+report['certificate_reactions']}
    steps = {s['id']:s for s in orientations(list(reactions.values()))}
    forbidden = set(report['forbidden_step_ids'])
    assert forbidden == set(parent['forbidden_step_ids']) | {
        r['id']+':hypothetical-right-to-left' for r in report['added_reactions']
        if r['hypothesis_type'] != 'glycerophospholipid-speciation'}
    evidence = json.loads(Path('data/reports/phase1-glycerophospholipid-synthesis.json').read_bytes())
    sources = {r['id']:r for r in evidence['reactions']}
    added = {r['id'] for r in report['added_reactions']}
    for r in report['added_reactions']:
        assert r['hypothesis_type'] in TYPES
        assert r == sources[r['id']]
    for r in reactions.values():
        assert balanced([r['left'],r['right']],cs)
    for cert in report['new_certificates']:
        assert not forbidden & {s['step_id'] for s in cert['steps']}
        assert added & {s['reaction_id'] for s in cert['steps']}
        validate_certificate(cert,steps,cs,set(report['external_exchange_compound_ids']),report['co2_compound_id'])
    old = {t['cannabisdb_id'] for t in parent['targets'] if t['net_status']=='exact-net-conversion-hypothesis'}
    covered = {t['cannabisdb_id'] for t in report['targets'] if t['net_status']=='exact-net-conversion-hypothesis'}
    assert old <= covered
    structures = {t['compound_id'] for t in report['targets'] if t['cannabisdb_id'] in covered-old}
    assert {c['compound_id'] for c in report['new_certificates']} == structures
    assert len(report['new_certificates']) == len(structures)
    assert report['summary']['new_certificate_records'] == len(covered-old)
    assert report['summary']['new_certificate_structures'] == len(structures)
    assert report['summary']['net_status_counts'] == dict(Counter(t['net_status'] for t in report['targets']))
    assert decode(encode(report,hashlib.sha256(path.read_bytes()).hexdigest())[::-1]) == report
