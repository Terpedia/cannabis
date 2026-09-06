import hashlib
import json
import re
from pathlib import Path
from rdkit import Chem
from cannabis_carbon.phase1_catalog import stable_id
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_reaction_completion_net import validate_certificate
from cannabis_carbon.phase1_odd_chain_supply import STARTER


def test_exact_carrier_matched_supply_panel_and_all_producers():
    root = Path('data/reports')
    read = lambda n: json.loads((root / ('phase1-' + n + '.json')).read_bytes())
    report = read('odd-chain-supply'); current = read('odd-chain-net')
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    layers = [json.loads((root / n).read_bytes()) for n in current['baseline_certificate_reports']]
    reactions, compounds, model = assemble_current(read('full-balanced-network'), read('marts-completions'), current, layers)
    steps = {s['id']: s for s in model.steps}
    assert report['external_exchange_compound_ids'] == current['external_exchange_compound_ids']
    assert report['forbidden_step_ids'] == current['forbidden_step_ids']
    assert [t['acyl_carbon_count'] for t in report['targets']] == list(range(3, 24, 2))
    tail = compounds[STARTER]['smiles'][17:]
    used = set()
    positive = 0
    for t in report['targets']:
        exact = Chem.MolToSmiles(Chem.MolFromSmiles('C' * t['acyl_carbon_count'] + tail), isomericSmiles=True)
        assert t['smiles'] == exact and t['compound_id'] == stable_id('structure', exact)
        cid = t['compound_id']
        assert t['present_in_model'] == (cid in compounds)
        producers = [s for s in model.steps if sum(p['coefficient'] for p in s['outputs'] if p['compound_id'] == cid)
                     > sum(p['coefficient'] for p in s['required_inputs'] if p['compound_id'] == cid)]
        assert t['producing_steps'] == producers
        used.update(s['reaction_id'] for s in producers)
        result = t['net_result']
        if result['status'] == 'exact-net-conversion-hypothesis':
            positive += 1
            validate_certificate({'compound_id': cid, **result}, steps, compounds,
                                 set(report['external_exchange_compound_ids']), report['co2_compound_id'])
            used.update(s['reaction_id'] for s in result['steps'])
        elif not producers:
            assert result['status'] == 'no-net-producing-candidate-equation'
        else:
            assert result['status'] in ('solver-reported-infeasible', 'solver-incomplete-or-failed',
                                        'numerical-solution-failed-exact-validation')
            assert 'solver_status' in result and 'solver_message' in result
    assert report['summary'] == {'assessed_exact_odd_chain_acyl_CoAs': 11,
                                'exact_net_certificates': positive, 'added_equations': 0,
                                'inventory_coverage_gain_claimed': 0}
    assert report['reactions'] == [reactions[r] for r in sorted(used)]
    participants = {p['compound_id'] for rid in used for side in ('left', 'right') for p in reactions[rid][side]}
    assert report['compounds'] == [compounds[c] for c in sorted(participants)]
