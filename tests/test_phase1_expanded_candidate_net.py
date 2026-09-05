import copy
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import pytest
from rdkit import Chem

from cannabis_carbon.phase1_expanded_candidate_net import build

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


@pytest.fixture(scope='module')
def inputs():
    return [read(n) for n in ('phase1-full-balanced-network', 'phase1-candidate-scope',
        'phase1-combined-catalog-evidence', 'phase1-candidate-net-flux', 'phase1-expanded-candidate-net')]


def test_candidate_union_inventory_and_unchanged_baseline(inputs):
    network, startup, supplement, baseline, report = inputs
    old = startup['candidate_reaction_evidence_ids']
    added = {e['reaction_id']: [e['id']] for e in supplement['enzyme_evidence']}
    assert not old.keys() & added.keys()
    assert report['candidate_reaction_evidence_ids'] == old | added
    assert len(report['candidate_reaction_evidence_ids']) == 1588
    assert [(t['cannabisdb_id'], t['compound_id']) for t in report['targets']] == [(t['cannabisdb_id'], t['compound_id']) for t in network['targets']]
    certificates = {c['compound_id']: c for c in report['certificates']}
    assert all(certificates[c['compound_id']] == c for c in baseline['certificates'])
    assert len(certificates) == 107
    assert len([t for t in report['targets'] if t['certificate_compound_id']]) == 108
    assert {t['cannabisdb_id'] for t in report['targets'] if t['new_net_certificate']} == {
        'CDB004791', 'CDB004808', 'CDB004818', 'CDB004839', 'CDB004887', 'CDB004953', 'CDB004992'}
    assert report['external_exchange_compound_ids'] == baseline['external_exchange_compound_ids']
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_every_equation_and_certificate_independently_replays(inputs):
    network, _, _, _, report = inputs
    original = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    reactions = {r['id']: r for r in report['reactions']}
    atoms = {}
    for rid, r in reactions.items():
        assert r['enzyme_evidence_ids'] == report['candidate_reaction_evidence_ids'][rid]
        balance = Counter()
        for side, sign in [('left', -1), ('right', 1)]:
            assert r[side] == original[rid][side]
            for p in r[side]:
                cid = p['compound_id']
                if cid not in atoms:
                    mol = Chem.AddHs(Chem.MolFromSmiles(compounds[cid]['smiles']))
                    atoms[cid] = Counter((a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms())
                    atoms[cid]['charge'] = Chem.GetFormalCharge(mol)
                for k, n in atoms[cid].items():
                    balance[k] += sign * Fraction(p['coefficient']) * n
        assert all(n == 0 for n in balance.values())
    exchange = set(report['external_exchange_compound_ids'])
    assert {compounds[c]['smiles'] for c in exchange if compounds[c]['carbon_count']} == {'O=C=O'}
    for cert in report['certificates']:
        net = defaultdict(Fraction)
        for step in cert['steps']:
            assert step['direction_mode'] in ('hypothetical-left-to-right', 'hypothetical-right-to-left')
            sign = 1 if step['direction_mode'] == 'hypothetical-left-to-right' else -1
            assert Fraction(step['extent']) > 0
            for side, multiplier in [('left', -sign), ('right', sign)]:
                for p in reactions[step['reaction_id']][side]:
                    net[p['compound_id']] += multiplier * Fraction(step['extent']) * Fraction(p['coefficient'])
        assert net[cert['compound_id']] == Fraction(cert['target_amount']) >= 1
        assert all(n >= 0 for c, n in net.items() if c not in exchange)
        assert cert['external_net_consumption'] == {c: str(-n) for c, n in net.items() if n < 0}
        assert cert['net_exports'] == {c: str(n) for c, n in net.items() if n > 0}
        assert cert['zero_net_internal_participants'] == sorted(c for c, n in net.items() if not n and c not in exchange)
        assert sum(n * compounds[c]['carbon_count'] for c, n in net.items()) == 0
        assert -net[report['co2_compound_id']] == Fraction(cert['net_carbon_in']) == Fraction(cert['net_carbon_out'])


def test_rebuild_nonmutation_and_duplicate_rejection(inputs, monkeypatch):
    network, startup, supplement, baseline, report = inputs
    # Solver paths can differ across supported SciPy/HiGHS versions. Replay
    # captured solver outcomes here; the test above independently proves every
    # certificate and full equation, without relying on optimizer reproducibility.
    outcomes = {t['compound_id']: {'status': t['net_status']} | {
        k: t[k] for k in ('solver_status', 'solver_message') if k in t} for t in report['targets']}
    for cert in report['certificates']:
        outcomes[cert['compound_id']] = {k: v for k, v in cert.items() if k not in ('compound_id', 'net_carbon_in', 'net_carbon_out')}
    monkeypatch.setattr('cannabis_carbon.phase1_expanded_candidate_net.NetModel.solve', lambda self, cid: outcomes[cid])
    before = copy.deepcopy((network, startup, supplement, baseline))
    rebuilt = build(network, startup, supplement, baseline)
    assert {k: v for k, v in rebuilt.items() if k != 'scipy_version'} == {k: v for k, v in report.items() if k not in ('source_sha256', 'scipy_version')}
    assert (network, startup, supplement, baseline) == before
    bad = copy.deepcopy(supplement)
    bad['enzyme_evidence'].append(bad['enzyme_evidence'][0])
    with pytest.raises(ValueError, match='Duplicate'):
        build(network, startup, bad, baseline)


def test_static_view_exact_sources_and_complete_evidence(inputs):
    report = inputs[-1]
    folder = ROOT / 'docs/data/expanded-net-view'
    manifest = json.loads((folder / 'index.json').read_text())
    payload = (folder / 'bundle.json').read_bytes()
    assert manifest['sha256'] == hashlib.sha256(payload).hexdigest()
    assert manifest['bytes'] == len(payload)
    for path, sha in manifest['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
    bundle = json.loads(payload)
    assert bundle['certificates'] == report['certificates']
    assert bundle['targets'] == report['targets']
    evidence = {e['id']: e for e in bundle['enzyme_evidence']}
    assert all(eid in evidence for r in bundle['reactions'] for eid in r['enzyme_evidence_ids'])
    assert len([r for r in bundle['reactions'] if r['is_new_catalog_candidate']]) == 5
