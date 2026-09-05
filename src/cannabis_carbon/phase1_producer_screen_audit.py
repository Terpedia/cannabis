"""Audit every exact producer gap against all versioned Phase 1 search rows."""
import copy
import hashlib
import json
from collections import Counter
from pathlib import Path

from .phase1_scope import write_rows


def build(audit, searches, core=None):
    targets = [t for t in audit['targets'] if t['exact_catalog_net_producer_reaction_ids']]
    ids = {rid for t in targets for rid in t['exact_catalog_net_producer_reaction_ids']}
    reactions = {r['id']: r for r in audit['reactions']}
    core_by_id = {r['id']: r for r in (core or {}).get('reactions', [])}
    rows = []
    for rid in sorted(ids):
        prior = [{'report': path, 'row': copy.deepcopy(r)} for path, report in sorted(searches.items())
                 for r in report['rows'] if r['reaction_id'] == rid]
        statuses = sorted({p['row']['search_status'] for p in prior})
        if not prior:
            disposition = 'no-prior-search-record'
        elif any(p['row'].get('passing_alignment_ids') for p in prior):
            disposition = 'prior-candidate-lead-not-admitted; review-evidence-context'
        elif 'weak-hits-only' in statuses:
            disposition = 'prior-weak-hits; improve-evidence-before-rescreen'
        elif 'no-hits' in statuses:
            disposition = 'prior-no-hits; not-biological-absence'
        elif statuses == ['no-reference-sequence']:
            disposition = 'reference-sequence-gap; inspect-discovery-and-retrieval'
        else:
            disposition = 'prior-status-needs-review'
        affected = [t for t in targets if rid in t['exact_catalog_net_producer_reaction_ids']]
        core_sources = sorted({s['source_reaction_id'] for s in reactions[rid].get('sources', [])} & core_by_id.keys())
        rows.append({'reaction_id': rid, 'reaction': copy.deepcopy(reactions[rid]),
            'target_ids': sorted(t['cannabisdb_id'] for t in affected),
            'identity_conflict_target_ids': sorted(t['cannabisdb_id'] for t in affected if t['source_identity_status'] == 'source-structure-disagreement'),
            'prior_searches': prior, 'prior_statuses': statuses, 'disposition': disposition,
            'core_source_reactions': [copy.deepcopy(core_by_id[s]) for s in core_sources],
            'core_evidence_boundary': 'Exact source-ID provenance links only. Legacy candidates, annotations and specialized searches are retained as recorded, not independently validated or admitted to the candidate model. A Phase 1 reference gap is not absence of legacy enzyme evidence.',
            'admitted_to_candidate_model': False})
    return {'schema': 'cannabis-carbon.phase1-producer-screen-audit.v1', 'rows': rows,
        'targets': copy.deepcopy(targets), 'search_report_paths': sorted(searches),
        'summary': {'target_records': len(targets), 'exact_catalog_producing_equations': len(rows),
            'prior_search_rows': sum(len(r['prior_searches']) for r in rows),
            'disposition_counts': dict(Counter(r['disposition'] for r in rows))},
        'claim_boundary': 'Exact catalog producer inventory, not complete CO2 routes. All original search rows and evidence classes are retained without promotion. A row with no reference sequence is not a completed negative proteome alignment; no hits and weak hits do not establish biological absence. Candidate hits from inferred equations require source-context review and are not silently admitted to the exact candidate model. Identity conflicts remain flagged. No new search, reaction, enzyme validation or atom tracing is claimed.'}


def run():
    source = Path('data/reports/phase1-no-producer-audit.json')
    paths = sorted(Path('data/reports').glob('phase1-*search.json'))
    searches = {str(p): json.loads(p.read_text()) for p in paths}
    if any('rows' not in r for r in searches.values()):
        raise ValueError('Unrecognized search report schema')
    core = Path('docs/data/networkdb.json')
    report = build(json.loads(source.read_text()), searches, json.loads(core.read_text()))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [source, core, *paths]}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-producer-screen-audit.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    records = [('metadata', 'report', {k: v for k, v in report.items() if k not in ('rows', 'targets')})]
    records.extend(('producer', r['reaction_id'], r) for r in report['rows'])
    records.extend(('target', r['cannabisdb_id'], r) for r in report['targets'])
    count = write_rows(records, sha, Path('data/derived/phase1-producer-screen-audit.ndjson'))
    print(json.dumps({'summary': report['summary'], 'rows': count, 'sha256': sha}))


if __name__ == '__main__':
    run()
