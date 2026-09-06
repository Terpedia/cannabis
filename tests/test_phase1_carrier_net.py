from cannabis_carbon.phase1_net_flux import NetModel


def step(name, left, right):
    return {'id': name, 'reaction_id': name, 'direction_mode': 'test-forward',
            'required_inputs': [{'compound_id': c, 'coefficient': 1} for c in left],
            'outputs': [{'compound_id': c, 'coefficient': 1} for c in right]}


def test_carrier_must_regenerate_without_import_or_accumulation():
    steps = [step('production', ['co2', 'red'], ['target', 'ox'])]
    model = NetModel([], {'co2'}, directed_steps=steps, conserved_ids={'red', 'ox'})
    assert model.solve('target')['status'] == 'solver-reported-infeasible'
    steps.append(step('regeneration', ['ox', 'light-reductant'], ['red']))
    model = NetModel([], {'co2', 'light-reductant'}, directed_steps=steps, conserved_ids={'red', 'ox'})
    result = model.solve('target')
    assert result['status'] == 'exact-net-conversion-hypothesis'
    assert result['external_net_consumption'] == {'co2': '1', 'light-reductant': '1'}
    assert set(result['zero_net_internal_participants']) == {'red', 'ox'}
    # No reverse orientations may be invented for explicit step input.
    assert [s['id'] for s in model.steps] == ['production', 'regeneration']


def test_full_network_keeps_all_explicit_steps_and_carrier_equalities():
    import json
    from pathlib import Path
    from cannabis_carbon.phase1_carrier_net import build
    from cannabis_carbon.phase1_net_flux import exact_net
    network = json.loads(Path('data/reports/phase1-carrier-network.json').read_bytes())
    model = build(network)
    assert model.steps == network['directed_steps']
    assert len(model.steps) == 33592
    assert len(model.conserved_ids) == 81
    assert model.conservation_matrix.shape == (81, 33592)
    assert not model.conserved_ids & model.exchange_ids
    matrix = model.matrix.tocsc()
    for j, s in enumerate(model.steps):
        expected = {model.index[c]: -float(v) for c, v in exact_net([s], ['1']).items()
                    if c in model.index and v}
        column = matrix.getcol(j)
        assert dict(zip(column.indices, column.data)) == expected
    assert (model.conservation_matrix != model.matrix[[model.index[c] for c in sorted(model.conserved_ids)]]).nnz == 0
