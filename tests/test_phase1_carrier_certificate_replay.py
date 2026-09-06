from cannabis_carbon.phase1_carrier_certificate_replay import replay


def step(name, consumed, product):
    return {'id': name, 'required_inputs': [{'compound_id': consumed, 'coefficient': 1}],
            'outputs': [{'compound_id': product, 'coefficient': 1}]}


def test_forced_depletion_is_not_hidden_by_carrier_choices():
    cert = {'compound_id': 'target', 'target_amount': '1', 'steps': [{'step_id': 'old', 'extent': '1'}]}
    r = replay(cert, {'old': [step('a', 'carrier', 'target')]}, set())
    assert r['unavoidable_internal_depletion'] == {'carrier': '1'}
    assert r['status'] == 'fixed-historical-flux-cannot-avoid-internal-depletion'
    r = replay(cert, {'old': [step('a', 'c1', 'target'), step('b', 'c2', 'target')]}, set())
    # Neither individual carrier is unavoidably depleted, but the alternatives
    # are not jointly feasible without a source. Do not label this a success.
    assert not r['unavoidable_internal_depletion']
    assert r['status'] == 'carrier-choice-allocation-unresolved'
    r = replay(cert, {'old': [step('a', 'co2', 'target')]}, {'co2'})
    assert r['status'] == 'fixed-flux-net-accounting-survives-not-full-pathway-validation'


def test_all_saved_certificates_and_inventory_are_audited():
    import hashlib
    import json
    from pathlib import Path
    from collections import Counter
    r = json.loads(Path('data/reports/phase1-carrier-certificate-replay.json').read_bytes())
    for path, sha in r['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    assert len(r['certificates']) == 2737
    assert len(r['targets']) == 6220
    assert len({c['compound_id'] for c in r['certificates']}) == 2737
    assert dict(Counter(c['status'] for c in r['certificates'])) == r['summary']['certificate_status_counts']
    assert r['summary']['qualified_co2_pathway_coverage'] is None


def test_every_reported_forced_depletion_has_exact_upper_bound():
    import json
    from pathlib import Path
    from fractions import Fraction
    from collections import defaultdict
    read = lambda n: json.loads(Path('data/reports/phase1-' + n + '.json').read_bytes())
    network = read('carrier-network'); report = read('carrier-certificate-replay')
    steps = {s['id']: s for s in network['directed_steps']}
    exchange = set(network['external_exchange_compound_ids'])
    for result in report['certificates']:
        upper = defaultdict(Fraction)
        for choice in result['step_choices']:
            alternatives = []
            for sid in choice['candidate_step_ids']:
                s = steps[sid]
                assert choice['original_step_id'] in s.get('inherited_allowed_step_ids', [s['id']])
                net = defaultdict(Fraction)
                for side, sign in (('required_inputs', -1), ('outputs', 1)):
                    for p in s[side]:
                        net[p['compound_id']] += sign * p['coefficient'] * Fraction(choice['extent'])
                alternatives.append(net)
            for cid in set().union(*(n.keys() for n in alternatives)):
                upper[cid] += max(n.get(cid, 0) for n in alternatives)
        forced = {cid: str(-v) for cid, v in upper.items() if v < 0 and cid not in exchange}
        assert result['unavoidable_internal_depletion'] == forced
