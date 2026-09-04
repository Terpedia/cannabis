import hashlib
import json
from pathlib import Path


def test_published_lazy_bundles_preserve_every_target_hypothesis_and_input():
    root = Path(__file__).resolve().parents[1]
    source = root / 'data/reports/phase1-target-hypotheses.json'
    report = json.loads(source.read_text())
    folder = root / 'docs/data/hypothesis-view'
    index = json.loads((folder / 'index.json').read_text())
    assert index['source_sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
    if index.get('enzyme_overlay_report'):
        from cannabis_carbon.phase1_screened_overlay import apply_overlay
        overlay = root / 'data/reports' / index['enzyme_overlay_report']
        assert index['enzyme_overlay_sha256'] == hashlib.sha256(overlay.read_bytes()).hexdigest()
        report = apply_overlay(report, json.loads(overlay.read_text()))
    assert index['summary'] == report['summary']
    assert [t['cannabisdb_id'] for t in index['targets']] == [t['cannabisdb_id'] for t in report['targets']]
    expected_hypotheses = {h['id']: h for h in report['hypotheses']}
    expected_reactions = {r['id']: r for r in report['reactions']}
    expected_compounds = {c['id']: c for c in report['compounds']}
    expected_evidence = {e['id']: e for e in report['enzyme_evidence']}
    loaded_targets, loaded_hypotheses = {}, {}
    for filename, metadata in index['files'].items():
        raw = (folder / filename).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == metadata['sha256']
        assert len(raw) == metadata['bytes']
        data = json.loads(raw)
        if filename.startswith('targets/'):
            loaded_targets[data['target']['cannabisdb_id']] = data
            assert [h['id'] for h in data['hypotheses']] == data['target']['hypothesis_ids']
            for h in data['hypotheses']:
                assert h['reaction_id'] == expected_hypotheses[h['id']]['reaction_id']
                assert h['has_candidate_enzyme_evidence'] == expected_hypotheses[h['id']]['has_candidate_enzyme_evidence']
        else:
            for rid, bundle in data.items():
                assert filename == f'reactions/{rid.split(":")[1][:2]}.json'
                assert bundle['reaction'] == expected_reactions[rid]
                for h in bundle['hypotheses']:
                    assert h == expected_hypotheses[h['id']]
                    assert h['id'] not in loaded_hypotheses
                    loaded_hypotheses[h['id']] = h
                cids = {m['compound_id'] for side in ('left', 'right') for m in bundle['reaction'][side]}
                assert {c['id'] for c in bundle['compounds']} == cids
                for c in bundle['compounds']:
                    assert {k: v for k, v in c.items() if k != 'labels'} == expected_compounds[c['id']]
                for e in bundle['enzyme_evidence']:
                    assert e == expected_evidence[e['id']]
    assert loaded_hypotheses == expected_hypotheses
    assert set(loaded_targets) == {t['cannabisdb_id'] for t in report['targets'] if t['hypothesis_ids']}
