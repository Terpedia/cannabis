import hashlib
import json
import pytest
from collections import Counter
from pathlib import Path
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_medium_boundary import UptakeLimitedModel, blocked_species
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_net_obstruction import validate as validate_obstruction


@pytest.mark.parametrize('name', ['recent-gains-medium', 'recent-gains-no-gases'])
def test_current_model_four_target_restriction_and_exact_proofs(name):
    root = Path('data/reports')
    read = lambda n: json.loads((root / ('phase1-' + n + '.json')).read_bytes())
    report = read(name); current = read('geranial-reduction-net')
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    layers = [json.loads((root / n).read_bytes()) for n in current['baseline_certificate_reports']]
    reactions, compounds, _ = assemble_current(read('full-balanced-network'), read('marts-completions'), current, layers)
    exchange = set(current['external_exchange_compound_ids'])
    blocked = blocked_species(compounds, exchange)
    if name == 'recent-gains-no-gases':
        blocked += [{'compound_id': c, 'smiles': compounds[c]['smiles'],
                     'reason': 'nitrogen-or-hydrogen-gas-uptake-diagnostic'}
                    for c in sorted(exchange) if compounds[c]['smiles'] in ('N#N', '[H][H]')]
    assert report['blocked_inputs'] == blocked
    uptake = exchange - {r['compound_id'] for r in blocked}
    if name == 'recent-gains-no-gases':
        assert not {compounds[c]['smiles'] for c in uptake} & {'N#N', '[H][H]'}
    assert set(report['allowed_uptake_compound_ids']) == uptake
    assert set(report['allowed_external_output_compound_ids']) == exchange
    assert report['forbidden_step_ids'] == current['forbidden_step_ids']
    assert {c for c in uptake if compounds[c]['carbon_count']} == {current['co2_compound_id']}
    model = UptakeLimitedModel(list(reactions.values()), exchange, uptake, current['forbidden_step_ids'])
    steps = {s['id']: s for s in model.steps}
    audited = {r['compound_id']: r for r in read('geranial-evidence-audit')['certificates']}
    assert len(report['targets']) == len(audited) == 4
    assert {r['compound_id'] for r in report['targets']} == audited.keys()
    used = set(exchange); used_reactions = set()
    for row in report['targets']:
        cid = row['compound_id']; prior = audited[cid]
        for key in ('cannabisdb_id', 'label', 'compound_id'):
            assert row[key] == prior[key]
        assert row['prior_certificate_blocked_inputs'] == sorted(set(prior['certificate']['external_net_consumption']) - uptake)
        result = row['restricted_net_result']; assert result['compound_id'] == cid
        if result['status'] == 'exact-net-conversion-hypothesis':
            validate_certificate(result, steps, compounds, uptake, current['co2_compound_id'])
            assert set(result['external_net_consumption']) <= uptake
            for s in result['steps']:
                used_reactions.add(s['reaction_id'])
                used.update(p['compound_id'] for side in ('required_inputs', 'outputs') for p in steps[s['step_id']][side])
        else:
            proof = row['obstruction_result']
            assert proof['compound_id'] == cid
            if proof['status'] == 'exact-stoichiometric-obstruction':
                assert validate_obstruction(model.model, cid, proof['weights']) == proof['weighted_steps']
                assert proof['checked_allowed_steps'] == len(steps)
                used.update(proof['weights'])
        used.add(cid)
    assert report['compounds'] == [compounds[c] for c in sorted(used)]
    assert report['certificate_reactions'] == [reactions[r] for r in sorted(used_reactions)]
    assert report['summary'] == {'assessed_target_records': 4, 'full_model_inventory_records': 6220,
        'balanced_equations': len(reactions), 'allowed_steps': len(steps),
        'net_status_counts': dict(Counter(r['restricted_net_result']['status'] for r in report['targets'])),
        'exact_obstructions': sum(r.get('obstruction_result', {}).get('status') == 'exact-stoichiometric-obstruction' for r in report['targets'])}
