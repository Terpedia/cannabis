"""Evidence-only overlay; completion chemistry and pathway metrics remain unchanged."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_completion_protein_discovery import BOUNDARY
from .phase1_new_protein_search import export_table as export_search


def build(completions, discovery, search):
    existing = {r['id']: r for r in discovery['existing_evidence_matches']}
    queued = {r['reaction_id']: r for r in discovery['rows']}
    searched = {r['reaction_id']: r for r in search['rows']}
    if len(searched) != len(search['rows']) or queued.keys() != searched.keys():
        raise ValueError('Search scope differs from discovery')
    alignments = {a['id']: a for a in search['passing_alignments']}
    rows = []
    for h in completions['completions']:
        eid, rid = h['id'], h['balanced_equation_id']
        if eid in existing:
            match = existing[eid]
            if match['reaction_id'] != rid or not match['evidence_ids']:
                raise ValueError('Invalid existing evidence match')
            rows.append({'id': eid, 'reaction_id': rid, 'category': 'existing-exact-equation-candidate-evidence',
                'existing_evidence_ids': match['evidence_ids'], 'screened_cannabis_proteins': [],
                'representative_alignments': [], 'has_candidate_lead': True,
                'claim_boundary': match['claim_boundary']})
            continue
        row = searched[rid]; original = queued[rid]
        if row['hypothesis_ids'] != [eid] or row['reference_matches'] != original['reference_matches']:
            raise ValueError('Search source reference join mismatch')
        allowed = {m['accession'] for m in row['reference_matches']}
        matches = [alignments[aid] for aid in row['passing_alignment_ids']]
        if any(a['reference_accession'] not in allowed or not a['passes_screen'] or
               a['identity_percent'] < 30 or a['query_coverage_percent'] < 50 or
               a['reference_coverage_percent'] < 50 or not 0 <= a['evalue'] <= 1e-5 for a in matches):
            raise ValueError('Alignment outside source references or thresholds')
        best = {}
        for a in sorted(matches, key=lambda a: (-a['bitscore'], a['id'])):
            best.setdefault(a['cannabis_accession'], a)
        if set(best) != set(row['screened_cannabis_proteins']):
            raise ValueError('Candidate protein set differs from passing alignments')
        rows.append({'id': eid, 'reaction_id': rid,
            'category': 'MARTS-source-homology-for-inferred-stoichiometry' if best else 'no-screened-candidate-lead',
            'has_candidate_lead': bool(best), 'existing_evidence_ids': [],
            'screened_cannabis_proteins': sorted(best), 'search_status': row['search_status'],
            'passing_alignment_ids': row['passing_alignment_ids'],
            'representative_alignments': [best[p] for p in sorted(best)],
            'representative_selection': 'Highest bit score per protein for display only, not a functional ranking; full search retains all passing alignments.',
            'missing_reference_sequences': row['reference_sequences_missing'],
            'validation_blockers': row['validation_blockers'], 'prior_source_reviews': original['prior_source_reviews'],
            'full_search_report': 'phase1-completion-protein-search.json', 'claim_boundary': BOUNDARY})
    if {r['id'] for r in rows} != {h['id'] for h in completions['completions']}:
        raise ValueError('Completion inventory changed')
    by_id = {r['id']: r for r in rows}
    return {'schema': 'cannabis-carbon.phase1-completion-protein-evidence.v1', 'rows': rows,
        'summary': {'completion_hypotheses': len(rows), 'category_counts': dict(Counter(r['category'] for r in rows)),
            'targets_with_any_candidate_lead': sum(any(by_id[hid]['has_candidate_lead'] for hid in t['completion_ids']) for t in completions['targets']),
            'targets_with_completions_but_no_candidate_lead': sum(bool(t['completion_ids']) and not any(by_id[hid]['has_candidate_lead'] for hid in t['completion_ids']) for t in completions['targets'])},
        'claim_boundary': 'Evidence sidecar only. Homology does not validate inferred stoichiometry, exact source-product identity, direction or a CO2 pathway. Existing candidate evidence is joined only on exact full equations. Atom tracing deferred.'}


def export(report_path, output, groups):
    raw = report_path.read_bytes(); report = json.loads(raw); digest = hashlib.sha256(raw).hexdigest()
    metadata = {k: v for k, v in report.items() if k not in {v for _, v, _ in groups}}
    count = 0
    with output.open('w') as handle:
        for kind, records, key in [('metadata', [metadata], None)] + [(k, report[v], key) for k, v, key in groups]:
            for row in records:
                handle.write(json.dumps({'record_kind': kind, 'record_id': row[key] if key else 'metadata',
                    'record_json': json.dumps(row, separators=(',', ':')), 'report_sha256': digest}) + '\n'); count += 1
    return count


def run():
    paths = [Path('data/reports', name + '.json') for name in ['phase1-marts-completions',
        'phase1-completion-protein-discovery', 'phase1-completion-protein-search']]
    reports = [json.loads(p.read_text()) for p in paths]
    if reports[2]['source_discovery_sha256'] != hashlib.sha256(paths[1].read_bytes()).hexdigest():
        raise ValueError('Search source digest mismatch')
    result = build(*reports)
    result['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    output = Path('data/reports/phase1-completion-protein-evidence.json')
    output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    counts = {}
    counts['discovery'] = export(paths[1], Path('data/derived/phase1-completion-protein-discovery.ndjson'),
        [('equation', 'rows', 'reaction_id'), ('existing_evidence', 'existing_evidence_matches', 'id'),
         ('excluded_reference', 'excluded_source_references', 'id')])
    counts['search'] = export_search(paths[2], Path('data/derived/phase1-completion-protein-search.ndjson'))
    counts['evidence'] = export(output, Path('data/derived/phase1-completion-protein-evidence.ndjson'), [('completion', 'rows', 'id')])
    for p in paths[1:] + [output]:
        Path('docs/data', p.name).write_bytes(p.read_bytes())
    print(json.dumps({'summary': result['summary'], 'export_rows': counts,
        'sha256': {p.stem: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths[1:] + [output]}}))


if __name__ == '__main__':
    run()
