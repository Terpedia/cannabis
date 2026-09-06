"""Inventory-wide impact of confirmed reactive-part/full-carrier identity conflicts."""
import hashlib
import json
from collections import Counter
from pathlib import Path


def run():
    root = Path('data/reports'); current_path = root / 'phase1-geranial-reduction-net.json'
    current = json.loads(current_path.read_bytes())
    paths = [current_path, Path('data/curation/rhea-carrier-identity-conflicts.json'),
             root / 'phase1-full-balanced-network.json'] + [root / n for n in current['baseline_certificate_reports']]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale source')
    conflicts = {r['reaction_id']: r for r in docs[1]['conflicts']}
    reactions = {r['id']: r for r in docs[2]['reactions']}
    compounds = {c['id']: c for c in current['compounds']}
    for rid, conflict in conflicts.items():
        r = reactions[rid]
        if set(conflict['source_reaction_ids']) != {s['source_reaction_id'] for s in r['sources']}:
            raise ValueError('Conflict does not match source reaction provenance')
        smiles = {compounds[p['compound_id']]['smiles'] for side in ('left', 'right') for p in r[side]}
        if not {c['incorrect_model_smiles'] for c in conflict['source_carriers']} <= smiles:
            raise ValueError('Expected reactive-part substitution absent')
    certs = {c['compound_id']: c for doc in [*docs[3:], current]
             for c in doc.get('certificates', []) + doc.get('new_certificates', [])}
    if set(certs) != {t['compound_id'] for t in current['targets'] if t['net_status'] == 'exact-net-conversion-hypothesis'}:
        raise ValueError('Incomplete certificate inventory')
    affected = {cid: [s for s in c['steps'] if s['reaction_id'] in conflicts] for cid, c in certs.items()}
    rows = []
    for t in current['targets']:
        cid = t['compound_id']; hits = affected.get(cid, [])
        rows.append({k: t[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'net_status')} | {
            'carrier_audit_status': 'confirmed-carrier-identity-conflict' if hits else
                'no-hit-in-two-reaction-screen-not-full-identity-validation' if cid in certs else 'no-saved-certificate-to-screen',
            'conflicting_steps': hits})
    counts = Counter(s['reaction_id'] for hits in affected.values() for s in hits)
    report = {'schema': 'cannabis-carbon.phase1-carrier-identity-audit.v1', 'targets': rows,
        'conflicts': list(conflicts.values()), 'model_reactions': [reactions[r] for r in conflicts],
        'summary': {'target_records': len(rows), 'screened_certificate_structures': len(certs),
            'affected_certificate_structures': sum(bool(s) for s in affected.values()),
            'affected_inventory_records': sum(bool(t['conflicting_steps']) for t in rows),
            'certificate_counts_by_conflicting_reaction': dict(counts), 'corrected_pathway_coverage_established': False},
        'claim_boundary': docs[1]['claim_boundary'],
        'remediation': 'Audit all Rhea composite/polymeric participants, retain full source carrier IDs and carrier-specific redox conservation, separate free-ion nutrient exchange, then recompute coverage. Do not invent full molecular formulas for unspecified proteins or merge distinct carrier classes. Source-faithful macromolecular equations require an explicit representation and balance policy. Quarantining identified projected equations is a sensitivity, not the completed correction.',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root / 'phase1-carrier-identity-audit.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
