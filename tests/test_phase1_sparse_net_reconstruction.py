from fractions import Fraction
import json
from pathlib import Path
import numpy as np
import pytest
from cannabis_carbon.phase1_net_flux import NetModel, exact_net, _reconstruct_sparse
from cannabis_carbon.phase1_scope import orientations


def reaction(rid, left, right):
    return {'id': rid,
            'left': [{'compound_id': c, 'coefficient': n} for c, n in left],
            'right': [{'compound_id': c, 'coefficient': n} for c, n in right]}


def test_sparse_reconstruction_matches_dense_with_fractional_and_zero_fluxes():
    reactions = [reaction('fix', [('CO2', 3), ('pool', 1)], [('A', 2)]),
                 reaction('release', [('A', 2)], [('pool', 1), ('target', 3)])]
    reactions += [reaction('unused-' + str(i), [('unused-A-' + str(i), 1)],
                            [('unused-B-' + str(i), 1)]) for i in range(1000)]
    steps = orientations(reactions)
    for active in (0.0, 1 / 3, 7 / 11, 1e-12):
        fluxes = np.zeros(len(steps))
        for i, s in enumerate(steps):
            if s['reaction_id'] in ('fix', 'release') and s['direction_mode'] == 'hypothetical-left-to-right':
                fluxes[i] = active
        amounts = [Fraction(float(x)).limit_denominator(1000000) for x in fluxes]
        dense = exact_net(steps, amounts)
        used, sparse = _reconstruct_sparse(steps, fluxes)
        assert used == [(s, n) for s, n in zip(steps, amounts) if n]
        assert {c: n for c, n in sparse.items() if n} == {c: n for c, n in dense.items() if n}
        participants = {m['compound_id'] for s, _ in used for side in ('required_inputs', 'outputs') for m in s[side]}
        assert {c: sparse[c] for c in participants} == {c: dense[c] for c in participants}


def test_sparse_validation_does_not_hide_bad_extents_or_source_coefficients():
    steps = orientations([reaction('r', [('A', 1)], [('B', 1)])])
    with pytest.raises(ValueError, match='Every step'):
        _reconstruct_sparse(steps, [])
    with pytest.raises(ValueError, match='Negative directed extent'):
        _reconstruct_sparse(steps, [-1, 0])
    # All reactions are validated at construction, even if later unused by a solve.
    with pytest.raises(ValueError, match='positive coefficients'):
        NetModel([reaction('bad', [('A', 0)], [('B', 1)])], {'CO2'})


def test_every_published_hydrolysis_certificate_keeps_exact_net_and_active_pools():
    report = json.loads(Path('data/reports/phase1-phosphatidate-hydrolysis-net.json').read_text())
    steps = orientations(report['certificate_reactions'])
    index = {s['id']: i for i, s in enumerate(steps)}
    for cert in report['new_certificates']:
        fluxes = np.zeros(len(steps))
        selected, extents = [], []
        for record in cert['steps']:
            i = index[record['step_id']]
            amount = Fraction(record['extent'])
            fluxes[i] = float(amount)
            selected.append(steps[i]); extents.append(amount)
        expected = exact_net(selected, extents)
        used, actual = _reconstruct_sparse(steps, fluxes)
        assert actual == expected
        assert len(used) == len(cert['steps'])
