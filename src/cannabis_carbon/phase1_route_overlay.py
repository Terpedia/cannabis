"""Join selected-route protein candidates into full-network and map evidence."""
import copy
import hashlib
import json
from pathlib import Path
from .phase1_screened_overlay import build_overlay, apply_overlay, export_table
from .phase1_scope import write_rows


def combine(parent, previous, added):
    allowed = {r['id'] for r in parent['reactions']}
    combined = copy.deepcopy(previous)
    if {e['id'] for e in previous['enzyme_evidence']} & {e['id'] for e in added['enzyme_evidence']}:
        raise ValueError('Evidence collision between screens')
    combined['enzyme_evidence'].extend(e for e in added['enzyme_evidence'] if e['reaction_id'] in allowed)
    combined['summary'] = {
        'equations_with_screened_evidence': len(combined['enzyme_evidence']),
        'distinct_cannabis_proteins': len({p['accession'] for e in combined['enzyme_evidence'] for p in e['screened_proteins']}),
        'route_screen_equations_in_map': sum(e['reaction_id'] in allowed for e in added['enzyme_evidence']),
        'route_screen_upstream_equations_outside_map': sum(e['reaction_id'] not in allowed for e in added['enzyme_evidence'])}
    combined['integrated_summary'] = apply_overlay(parent, combined)['summary']
    return combined


def annotate_routes(certificates, added):
    evidence = {e['reaction_id']: e for e in added['enzyme_evidence']}
    routes = []
    for route in certificates['routes']:
        steps = []
        for step in route['steps']:
            new = evidence.get(step['reaction_id'])
            ids = step['enzyme_evidence_ids'] + ([new['id']] if new else [])
            steps.append({'step_id': step['id'], 'reaction_id': step['reaction_id'],
                'enzyme_evidence_ids': ids, 'has_candidate_enzyme_evidence': bool(ids),
                'new_screened_protein_ids': [p['accession'] for p in new['screened_proteins']] if new else []})
        missing = [s['step_id'] for s in steps if not s['has_candidate_enzyme_evidence']]
        routes.append({'cannabisdb_id': route['cannabisdb_id'], 'scenario_id': route['scenario_id'],
            'compound_id': route['compound_id'], 'steps': steps, 'missing_enzyme_step_ids': missing,
            'first_missing_enzyme_step_id': missing[0] if missing else None,
            'status': route['status'], 'biological_status': route['biological_status']})
    return {'schema': 'cannabis-carbon.phase1-route-evidence-status.v1', 'routes': routes,
        'summary': {'routes': len(routes), 'routes_with_missing_candidate_evidence': sum(bool(r['missing_enzyme_step_ids']) for r in routes),
            'distinct_equations_without_candidate_evidence': len({s['reaction_id'] for r in routes for s in r['steps'] if not s['has_candidate_enzyme_evidence']})},
        'claim_boundary': 'Evidence-only overlay. Original certificate equations, extents, seeds, direction assumptions and biological status are unchanged. Missing candidate evidence is not biological absence; candidates are not confirmed activity.'}


def run():
    names = ['phase1-target-hypotheses', 'phase1-screened-enzyme-overlay', 'phase1-full-balanced-network',
             'phase1-route-certificates', 'phase1-route-protein-search', 'phase1-route-references']
    paths = {name: Path('data/reports', name + '.json') for name in names}
    inputs = {name: json.loads(p.read_text()) for name, p in paths.items()}
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()}
    parent, previous, network, cert, search, discovery = [inputs[n] for n in names]
    checks = [previous['source_sha256'], cert['source_sha256'], discovery['source_sha256']]
    for check in checks:
        for path, digest in check.items():
            if path in hashes and hashes[path] != digest:
                raise ValueError('Overlay input checksum mismatch')
    if search['source_discovery_sha256'] != hashes[str(paths['phase1-route-references'])]:
        raise ValueError('Search discovery checksum mismatch')
    added = build_overlay({'reactions': network['reactions'], 'hypotheses': []}, search, 'phase1-route-protein-search.json')
    if any(r['enzyme_evidence_ids'] for r in network['reactions'] if r['id'] in {e['reaction_id'] for e in added['enzyme_evidence']}):
        raise ValueError('Route screen unexpectedly overlaps prior evidence')
    combined = combine(parent, previous, added)
    status = annotate_routes(cert, added)
    for name, report in [('phase1-route-enzyme-overlay', added), ('phase1-combined-enzyme-overlay', combined), ('phase1-route-evidence-status', status)]:
        report['source_sha256'] = hashes
        payload = json.dumps(report, separators=(',', ':')) + '\n'
        for folder in ('data/reports', 'docs/data'):
            Path(folder, name + '.json').write_text(payload)
        source = Path('data/reports', name + '.json')
        output = source.with_suffix('.ndjson')
        if 'enzyme_evidence' in report:
            export_table(source, output)
        else:
            rows = [('route', r['scenario_id'] + ':' + r['cannabisdb_id'], r) for r in report['routes']]
            rows.append(('metadata', 'status', {k: v for k, v in report.items() if k != 'routes'}))
            write_rows(rows, hashlib.sha256(payload.encode()).hexdigest(), output)
        print(name, json.dumps(report.get('integrated_summary', report['summary'])), flush=True)


if __name__ == '__main__':
    run()
