import hashlib
import json
from pathlib import Path

from cannabis_carbon.phase1_current_reaction_gaps import build, producing_steps
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current


def test_net_producers_sum_repeated_participants_and_obey_direction_exclusions():
    reaction = {'id': 'test:one', 'left': [
        {'compound_id': 'a', 'coefficient': 1}, {'compound_id': 'a', 'coefficient': 1},
        {'compound_id': 'b', 'coefficient': 1}], 'right': [
        {'compound_id': 'a', 'coefficient': 2}, {'compound_id': 'c', 'coefficient': 1}]}
    unrestricted = producing_steps([reaction])
    assert not unrestricted['a']
    assert unrestricted['b'] == {'test:one:hypothetical-right-to-left'}
    assert unrestricted['c'] == {'test:one:hypothetical-left-to-right'}
    restricted = producing_steps([reaction], ['test:one:hypothetical-right-to-left'])
    assert not restricted['a'] and not restricted['b']
    assert restricted['c'] == unrestricted['c']


def test_every_current_gap_and_diagnostic_is_reproducible_from_full_model():
    root = Path('data/reports')
    read = lambda p: json.loads(p.read_bytes())
    report = read(root/'phase1-current-reaction-gaps.json')
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    current = read(root/'phase1-glycerophospholipid-net.json')
    reactions, compounds, _ = assemble_current(
        read(root/'phase1-full-balanced-network.json'), read(root/'phase1-marts-completions.json'),
        current, [read(root/n) for n in current['baseline_certificate_reports']])
    before = json.dumps(current, sort_keys=True)
    rebuilt = build(current, compounds, reactions, read(root/'phase1-no-producer-audit.json'))
    assert rebuilt == {k:v for k,v in report.items() if k != 'source_sha256'}
    assert json.dumps(current, sort_keys=True) == before
    expected = {t['cannabisdb_id']:t for t in current['targets']
                if t['net_status'] == 'no-net-producing-equation'}
    assert len(report['targets']) == len(expected) == 3097
    assert report['summary']['balanced_equations'] == 18082
    assert report['summary']['coverage_gain_claimed'] == 0
    assert sum(report['summary']['category_counts'].values()) == len(expected)
    for row in report['targets']:
        original = expected[row['cannabisdb_id']]
        for key in ('compound_id', 'net_status', 'label'):
            assert row[key] == original[key]
        assert row['canonical_smiles'] == compounds[row['compound_id']]['smiles']
        assert all(row['compound_id'] not in values for values in row['diagnostic_alternatives'].values())
        assert set(row['excluded_producing_step_ids']) <= set(current['forbidden_step_ids'])
