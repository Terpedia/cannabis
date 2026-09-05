import copy
import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_balance_reference import concrete_participants
from cannabis_carbon.phase1_catalog import stable_id
from cannabis_carbon.phase1_marts_completions import build, balanced, composition
from cannabis_carbon.phase1_completion_view import build as view_build
from cannabis_carbon.balance import _reaction_smiles_balance


def network(smiles):
    compounds, sides = {}, []
    for side in concrete_participants(smiles):
        members = []
        for p in side:
            cid = stable_id('structure', p['smiles'])
            compounds[cid] = {'id': cid, **{k: v for k, v in p.items() if k != 'coefficient'}}
            members.append({'compound_id': cid, 'coefficient': p['coefficient']})
        sides.append(members)
    return {'compounds': list(compounds.values()), 'reactions': [{'id': 'R', 'left': sides[0], 'right': sides[1],
        'sources': [{'source_layer': 'terpedia-full-rhea-catalog'}]}]}


def audit(smiles):
    return {'source_ledger': [{'id': 'S', 'balance_status': 'imbalanced', 'source_record': {'reaction_smarts': smiles}}], 'targets': []}


def test_formula_template_preserves_original_product_and_flags_inference():
    result = build(audit('C=C>>COC'), network('C=C.O>>CCO'))
    assert len(result['completions']) == 1
    h = result['completions'][0]
    cs = {c['id']: c for c in result['compounds']}
    assert cs[h['original_organic_pair'][1]]['smiles'] == 'COC'
    assert h['reference_templates'][0]['product_match'].startswith('composition-only')
    assert h['enzyme_evidence_ids'] == []
    assert balanced([h['left'], h['right']], cs)
    added = h['inferred_inorganic_participants_in_MARTS_orientation']
    assert [cs[p['compound_id']]['smiles'] for p in added[0]] == ['O']
    assert added[1] == []


def test_stereo_isotope_charge_and_all_inputs_are_not_relaxed():
    assert not build(audit('C[C@@H](Cl)Br>>C[C@@H](Cl)Br'), network('C[C@H](Cl)Br.O>>C[C@H](Cl)Br.O'))['completions']
    assert not build(audit('C=C>>[13CH3]OC'), network('C=C.O>>CCO'))['completions']
    assert not build(audit('C=C>>CC[O-]'), network('C=C.O>>CCO'))['completions']
    assert not build(audit('C=C.C>>CCO'), network('C=C.O>>CCO'))['completions']
    assert not build(audit('C=C>>CCO'), network('C=C.CO>>CCO.C'))['completions']
    assert composition('[13CH3]O') != composition('CO')


def test_published_completions_preserve_all_sources_and_replay_templates():
    root = Path(__file__).resolve().parents[1]
    raw = (root / 'data/reports/phase1-marts-completions.json').read_bytes()
    report = json.loads(raw)
    assert raw == (root / 'docs/data/phase1-marts-completions.json').read_bytes()
    inputs = []
    for name, digest in report['source_sha256'].items():
        data = (root / name).read_bytes(); assert hashlib.sha256(data).hexdigest() == digest
        inputs.append(json.loads(data))
    before = copy.deepcopy(inputs)
    assert build(*inputs) == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert inputs == before
    assert len(report['targets']) == 6220
    assert {sid for v in report['variants'] for sid in v['source_record_ids']} == {s['id'] for s in inputs[0]['source_ledger'] if s['balance_status'] == 'imbalanced'}
    cs = {c['id']: c for c in report['compounds']}
    refs = {r['id']: r for r in report['reference_reactions']}
    for h in report['completions']:
        assert balanced([h['left'], h['right']], cs)
        smiles = '>>'.join('.'.join(cs[p['compound_id']]['smiles'] for p in h[side] for _ in range(p['coefficient'])) for side in ('left', 'right'))
        element, charge = _reaction_smiles_balance(smiles)
        assert element['status'] == charge['status'] == 'balanced'
        assert h['enzyme_evidence_ids'] == []
        forward = h['marts_forward_direction'] == 'hypothetical-left-to-right'
        sides = [h['left'], h['right']] if forward else [h['right'], h['left']]
        for i, side in enumerate(sides):
            organic = [p for p in side if cs[p['compound_id']]['carbon_count']]
            assert organic == [{'compound_id': h['original_organic_pair'][i], 'coefficient': 1}]
        for t in h['reference_templates']:
            ref = refs[t['reference_reaction_id']]
            rs = [ref['left'], ref['right']] if t['reference_direction'] == 'hypothetical-left-to-right' else [ref['right'], ref['left']]
            assert [[p for p in side if not cs[p['compound_id']]['carbon_count']] for side in rs] == h['inferred_inorganic_participants_in_MARTS_orientation']


def test_compact_view_has_exact_source_manifest_and_payload():
    root = Path(__file__).resolve().parents[1]; folder = root / 'docs/data/completion-view'
    manifest = json.loads((folder / 'index.json').read_text()); raw = (folder / manifest['file']).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == manifest['sha256'] and len(raw) == manifest['bytes']
    sources = []
    for name, digest in manifest['source_sha256'].items():
        data = (root / name).read_bytes(); assert hashlib.sha256(data).hexdigest() == digest
        sources.append(json.loads(data))
    assert json.loads(raw) == view_build(*sources)
    assert len(raw) < 6_000_000
