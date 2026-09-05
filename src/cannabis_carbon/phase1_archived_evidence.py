"""Compact archived-source evidence supplement; original reports remain immutable."""
import copy
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_scope import write_rows


def build(discovery, search, parent):
    sources = {r['reaction_id']: r for r in discovery['rows']}
    prior = {r['reaction_id']: r for r in parent['rows']}
    if len(search['rows']) != len(sources) or {r['reaction_id'] for r in search['rows']} != sources.keys():
        raise ValueError('Archived search scope mismatch')
    alignments = {a['id']: a for a in search['passing_alignments']}
    rows = []
    for r in search['rows']:
        source = sources[r['reaction_id']]
        if r['reference_matches'] != source['reference_matches'] or r['hypothesis_ids'] != source['hypothesis_ids']:
            raise ValueError('Archived reference join mismatch')
        old = prior[r['reaction_id']]
        if r['hypothesis_ids'] != [old['id']]:
            raise ValueError('Archived completion identity mismatch')
        allowed = {m['accession'] for m in source['reference_matches']}
        matches = [alignments[a] for a in r['passing_alignment_ids']]
        if any(a['reference_accession'] not in allowed or not a['passes_screen'] or
               a['identity_percent'] < 30 or a['query_coverage_percent'] < 50 or
               a['reference_coverage_percent'] < 50 or not 0 <= a['evalue'] <= 1e-5 for a in matches):
            raise ValueError('Invalid archived alignment')
        best = {}
        for a in sorted(matches, key=lambda a: (-a['bitscore'], a['id'])):
            best.setdefault(a['cannabis_accession'], a)
        if set(best) != set(r['screened_cannabis_proteins']):
            raise ValueError('Archived candidate set mismatch')
        rows.append({'id': old['id'], 'reaction_id': r['reaction_id'],
            'has_archived_candidate_lead': bool(best), 'previous_candidate_lead': old['has_candidate_lead'],
            'target_ids': r['target_ids'], 'search_status': r['search_status'],
            'screened_cannabis_proteins': sorted(best), 'passing_alignment_ids': r['passing_alignment_ids'],
            'representative_alignments': [best[p] for p in sorted(best)],
            'archive_resolution_links': source['archive_resolution_links'],
            'prior_source_reviews': source['prior_source_reviews'],
            'validation_blockers': r['validation_blockers'], 'evidence_class': r['evidence_class'],
            'proposed_test': r['proposed_test'], 'source_report': 'phase1-archived-protein-search.json',
            'claim_boundary': discovery['claim_boundary']})
    return {'schema': 'cannabis-carbon.phase1-archived-evidence.v1', 'rows': rows,
        'summary': {'equations_reviewed': len(rows), 'equations_with_archive_candidates': sum(r['has_archived_candidate_lead'] for r in rows),
            'new_equations_with_candidate_lead': sum(r['has_archived_candidate_lead'] and not r['previous_candidate_lead'] for r in rows),
            'search_status_counts': dict(Counter(r['search_status'] for r in rows))},
        'claim_boundary': 'Supplemental original-source sequence evidence only; not confirmed activities or a new CO2 pathway. Prior completion chemistry, candidate-scope, and net-connectivity reports are unchanged. Atom tracing remains deferred.'}


def apply(parent, supplement):
    result = copy.deepcopy(parent)
    rows = {r['id']: r for r in result['rows']}
    for extra in supplement['rows']:
        row = rows[extra['id']]
        if row['reaction_id'] != extra['reaction_id'] or row['has_candidate_lead'] != extra['previous_candidate_lead']:
            raise ValueError('Archive supplement parent mismatch')
        row['archived_source_screen'] = extra
        if not extra['has_archived_candidate_lead']:
            continue
        row['prior_evidence_category'] = row['category']
        row['category'] = 'combined-original-MARTS-source-homology' if row['has_candidate_lead'] else 'original-MARTS-archived-source-homology'
        row['has_candidate_lead'] = True
        row['screened_cannabis_proteins'] = sorted(set(row['screened_cannabis_proteins']) | set(extra['screened_cannabis_proteins']))
        row['passing_alignment_ids'] = sorted(set(row.get('passing_alignment_ids', [])) | set(extra['passing_alignment_ids']))
        best = {}
        for a in sorted(row['representative_alignments'] + extra['representative_alignments'], key=lambda a: (-a['bitscore'], a['id'])):
            best.setdefault(a['cannabis_accession'], a)
        row['representative_alignments'] = [best[p] for p in sorted(best)]
        row['validation_blockers'] = sorted(set(row.get('validation_blockers', [])) | set(extra['validation_blockers']))
    result['summary']['category_counts'] = dict(Counter(r['category'] for r in result['rows']))
    # Target coverage is recomputed by the viewer from exact completion membership.
    result['summary']['archive_supplement'] = supplement['summary']
    return result


def run():
    paths = [Path('data/reports', n + '.json') for n in ['phase1-archived-references',
        'phase1-archived-protein-search', 'phase1-completion-protein-evidence']]
    inputs = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    if inputs[1]['source_discovery_sha256'] != hashes[str(paths[0])]:
        raise ValueError('Archive search discovery checksum mismatch')
    report = build(*inputs); report['source_sha256'] = hashes
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    output = Path('data/reports/phase1-archived-evidence.json'); output.write_text(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k != 'rows'})]
    rows.extend(('completion', r['id'], r) for r in report['rows'])
    count = write_rows(rows, digest, Path('data/derived/phase1-archived-evidence.ndjson'))
    for p in paths[:2] + [output]:
        Path('docs/data', p.name).write_bytes(p.read_bytes())
    print(json.dumps({'summary': report['summary'], 'rows': count, 'sha256': digest}))


if __name__ == '__main__':
    run()
