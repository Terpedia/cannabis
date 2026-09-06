import hashlib
import json
from pathlib import Path


def test_complete_inventory_carrier_conflict_impact_without_overclaim():
    root = Path('data/reports')
    read = lambda n: json.loads((root / ('phase1-' + n + '.json')).read_bytes())
    report = read('carrier-identity-audit'); current = read('geranial-reduction-net')
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    evidence = json.loads(Path('data/curation/rhea-carrier-identity-conflicts.json').read_bytes())
    assert report['conflicts'] == evidence['conflicts']
    assert len({c['id'] for r in report['conflicts'] for c in r['source_carriers']}) == 4
    assert all(c['id'].startswith('RHEA-COMP:') for r in report['conflicts'] for c in r['source_carriers'])
    conflicts = {c['reaction_id'] for c in report['conflicts']}
    network = {r['id']: r for r in read('full-balanced-network')['reactions']}
    assert report['model_reactions'] == [network[c['reaction_id']] for c in report['conflicts']]
    layers = [json.loads((root / n).read_bytes()) for n in current['baseline_certificate_reports']] + [current]
    certs = {c['compound_id']: c for r in layers for c in r.get('certificates', []) + r.get('new_certificates', [])}
    assert len(certs) == report['summary']['screened_certificate_structures'] == 2737
    assert len(report['targets']) == len(current['targets']) == 6220
    affected = set(); counts = {rid: 0 for rid in conflicts}
    for row, target in zip(report['targets'], current['targets']):
        for key in ('cannabisdb_id', 'label', 'compound_id', 'net_status'):
            assert row[key] == target[key]
        cid = row['compound_id']; hits = [s for s in certs.get(cid, {}).get('steps', []) if s['reaction_id'] in conflicts]
        assert row['conflicting_steps'] == hits
        if hits:
            assert row['carrier_audit_status'] == 'confirmed-carrier-identity-conflict'
            affected.add(cid)
        elif cid in certs:
            assert row['carrier_audit_status'] == 'no-hit-in-two-reaction-screen-not-full-identity-validation'
        else:
            assert row['carrier_audit_status'] == 'no-saved-certificate-to-screen'
    for cert in certs.values():
        for s in cert['steps']:
            if s['reaction_id'] in conflicts:
                counts[s['reaction_id']] += 1
    assert report['summary']['certificate_counts_by_conflicting_reaction'] == counts
    assert len(affected) == report['summary']['affected_certificate_structures'] == 2479
    assert sum(bool(r['conflicting_steps']) for r in report['targets']) == report['summary']['affected_inventory_records'] == 2481
    assert not report['summary']['corrected_pathway_coverage_established']
