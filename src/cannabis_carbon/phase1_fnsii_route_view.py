"""Static Cytoscape view of the four preserved FNS-II sensitivity scenarios."""
import hashlib
import json
from pathlib import Path
from .phase1_net_view import build as attach_evidence

SOURCES = ('data/reports/phase1-fnsii-addition-sensitivity.json',
    'docs/data/remaining-net-view/bundle.json',
    'data/reports/phase1-fnsii-redox-hypothesis.json',
    'data/reports/phase1-target-hypotheses.json',
    'data/reports/phase1-screened-enzyme-overlay.json',
    'data/reports/phase1-route-enzyme-overlay.json',
    'data/reports/phase1-expanded-candidate-net.json',
    'data/reports/phase1-purine-candidate-net.json',
    'data/reports/phase1-replacement-candidate-net.json',
    'data/reports/phase1-thiolase-candidate-net.json',
    'data/reports/phase1-remaining-candidate-net.json')


def build():
    report, existing, redox, *evidence_sources = [json.loads(Path(p).read_text()) for p in SOURCES]
    for document in (report, redox):
        for path, sha in document['source_sha256'].items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed sensitivity source')
    scenarios = []
    for s in report['scenarios']:
        certificates, targets = [], []
        for row in s['rows']:
            feasible = row['status'] == 'exact-net-conversion-hypothesis'
            if feasible:
                certificates.append(row)
            targets.append({k: row[k] for k in ('cannabisdb_id', 'compound_id', 'label')} |
                {'net_status': row['status'], 'startup_status': 'unestablished; pre-existing internal pools may be required',
                 'certificate_compound_id': row['compound_id'] if feasible else None,
                 'missing_candidate_reaction_ids': s['added_reaction_ids'] if feasible else []})
        scenarios.append({**{k: s[k] for k in ('id', 'forbidden_step_ids', 'added_reaction_ids', 'reaction_count')},
            'targets': targets, 'certificates': certificates,
            'summary': {'target_records': len(targets), 'target_status_counts': {
                'exact-net-conversion-hypothesis': len(certificates), **s['summary']}}})
    used = {s['reaction_id'] for scenario in scenarios for c in scenario['certificates'] for s in c['steps']}
    additions = {report['hypothetical_chi_reaction_id'], report['hypothetical_fnsii_reaction_id']}
    reactions = []
    for r in report['reactions']:
        if r['id'] not in used:
            continue
        is_added = r['id'] in additions
        ids = report['baseline_candidate_reaction_evidence_ids'].get(r['id'], [])
        if not is_added and not ids:
            raise ValueError('Missing baseline candidate evidence')
        item = {**r, 'enzyme_evidence_ids': ids, 'missing_candidate_evidence': is_added,
            'is_route_sensitivity': is_added}
        if is_added:
            item['hypothesis_assumptions'] = (
                redox['unverified_requirements'] if r['id'] == report['hypothetical_fnsii_reaction_id']
                else ['Exact Cannabis CHI activity and stereoselectivity remain unverified'])
            item['sources'] = [*r['sources'], {'source_urls': [
                'https://terpedia.github.io/cannabis/data/fnsii-route-sensitivity-bundle.json']}]
        reactions.append(item)
    required = {p['compound_id'] for r in reactions for side in ('left', 'right') for p in r[side]}
    required.update(c for s in scenarios for cert in s['certificates']
                    for c in (*cert['external_net_consumption'], *cert['net_exports']))
    base = {'schema': 'cannabis-fnsii-route-view-v1', 'view_scenario': 'fnsii-route-sensitivity',
        'model_eligible': False, 'scenario_options': scenarios,
        **next(s for s in scenarios if s['id'] == 'CHI-and-FNSII'),
        'compounds': [c for c in report['compounds'] if c['id'] in required], 'reactions': reactions,
        'view_boundary': 'CHI/FNS-II sensitivity for two selected flavonoids. Red edges mark the assumed enzyme steps. Both additions are needed for conditional net conversion from CO2. Carrier partnership and native activity remain unverified.',
        'claim_boundary': report['claim_boundary'],
        'source_sha256': {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in SOURCES}}
    bundle = attach_evidence(base, [existing, *evidence_sources])
    # Reuse established labels wherever available without changing identities.
    labels = {c['id']: c.get('labels', []) for c in existing['compounds']}
    for c in bundle['compounds']:
        if not c['labels']:
            c['labels'] = labels.get(c['id'], [])
    return bundle


if __name__ == '__main__':
    bundle = build()
    payload = json.dumps(bundle, separators=(',', ':')) + '\n'
    folder = Path('docs/data/fnsii-net-view')
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'bundle.json').write_text(payload)
    manifest = {'schema': bundle['schema'], 'file': 'bundle.json',
        'sha256': hashlib.sha256(payload.encode()).hexdigest(), 'bytes': len(payload.encode()),
        'source_sha256': bundle['source_sha256']}
    (folder / 'index.json').write_text(json.dumps(manifest, separators=(',', ':')) + '\n')
    print(json.dumps({'bytes': manifest['bytes'], 'reactions': len(bundle['reactions']),
        'compounds': len(bundle['compounds']), 'evidence': len(bundle['enzyme_evidence'])}))
