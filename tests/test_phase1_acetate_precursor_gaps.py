import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current
from cannabis_carbon.phase1_net_obstruction import validate


def test_four_precursor_gaps_and_exact_target_obstructions():
    root = Path('data/reports')
    read = lambda n: json.loads((root / ('phase1-' + n + '.json')).read_bytes())
    report = read('acetate-precursor-gaps'); current = read('alcohol-acetates-net')
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    layers = [json.loads((root / n).read_bytes()) for n in current['baseline_certificate_reports']]
    reactions, compounds, model = assemble_current(read('full-balanced-network'), read('marts-completions'), current, layers)
    assert {t['cannabisdb_id'] for t in report['targets']} == {'CDB000192','CDB000283','CDB000321','CDB000585'}
    assert report['external_exchange_compound_ids'] == current['external_exchange_compound_ids']
    assert report['forbidden_step_ids'] == current['forbidden_step_ids']
    assert report['summary'] == {'assessed_blocked_records': 4, 'exact_target_obstructions': 4,
        'alcohols_without_producing_steps': 1, 'balanced_equations': len(reactions), 'allowed_steps': len(model.steps), 'coverage_change': 0}
    used = set()
    steps = {s['id']: s for s in model.steps}
    for row in report['targets']:
        cid = row['alcohol_compound_id']
        assert row['alcohol_smiles'] == compounds[cid]['smiles']
        expected = [s for s in model.steps if sum(p['coefficient'] for p in s['outputs'] if p['compound_id'] == cid) >
                    sum(p['coefficient'] for p in s['required_inputs'] if p['compound_id'] == cid)]
        assert expected == row['alcohol_producing_steps']
        assert row['alcohol_net_result']['status'] == ('solver-reported-infeasible' if expected else 'no-net-producing-candidate-equation')
        proof = row['target_obstruction']
        assert proof['status'] == 'exact-stoichiometric-obstruction'
        assert proof['weighted_steps'] == validate(model,row['compound_id'],proof['weights'])
        assert proof['checked_allowed_steps'] == len(model.steps)
        assert cid in proof['weights']
        used.add(row['acetylation_hypothesis_id'])
        used.update(s['reaction_id'] for s in expected)
        used.update(steps[s['step_id']]['reaction_id'] for s in proof['weighted_steps'])
    assert report['reactions'] == [reactions[r] for r in sorted(used)]
    cids = {p['compound_id'] for rid in used for side in ('left','right') for p in reactions[rid][side]}
    cids.update(c for row in report['targets'] for c in row['target_obstruction']['weights'])
    assert report['compounds'] == [compounds[c] for c in sorted(cids)]
