"""Prior-screen audit and reviewed references for alternative purine net gaps."""
import copy
import hashlib
import json
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .phase1_new_references import attach
from .phase1_reference_discovery import direction_families, exact_annotations
from .phase1_scope import write_rows

PRIOR = ('phase1-new-protein-search', 'phase1-route-protein-search', 'phase1-catalog-protein-search',
         'phase1-backfill-protein-search', 'phase1-completion-protein-search', 'phase1-archived-protein-search',
         'phase1-targeted-protein-search', 'phase1-family-protein-search', 'phase1-plant-purine-search')


def queue(report, prior, families):
    if len({g['reaction_id'] for g in report['catalog_candidate_gaps']}) != len(report['catalog_candidate_gaps']):
        raise ValueError('Duplicate candidate gap')
    audit, rows = [], []
    for gap in report['catalog_candidate_gaps']:
        rid = gap['reaction_id']
        previous = [{'report': name, 'search_status': r.get('search_status'),
                     'reference_sequences_present': r.get('reference_sequences_present', []),
                     'passing_alignment_ids': r.get('passing_alignment_ids', [])}
                    for name, p in prior.items() for r in p.get('rows', []) if r.get('reaction_id') == rid]
        if any(r['passing_alignment_ids'] for r in previous):
            raise ValueError('Passing prior candidate absent from model')
        audit.append({**copy.deepcopy(gap), 'prior_screens': previous,
                      'disposition': 'retain-prior-search-result' if previous else 'new-reference-discovery'})
        if previous:
            continue
        sources = sorted({s['source_reaction_id'].upper() for s in gap['source_joins']})
        probes = {u['probe_compound_id'] for u in gap['selected_uses']}
        targets = sorted({t['cannabisdb_id'] for t in report['focused_targets'] if t['compound_id'] in probes})
        rows.append({'reaction_id': rid, 'left': copy.deepcopy(gap['left']), 'right': copy.deepcopy(gap['right']),
            'sources': copy.deepcopy(gap['source_joins']), 'source_reaction_ids': sources,
            'selected_probe_compound_ids': sorted(probes), 'selected_uses': copy.deepcopy(gap['selected_uses']),
            'target_ids': targets, 'priority_target_ids': targets, 'hypothesis_ids': [],
            'rhea_families': {sid: families[sid] for sid in sources if sid in families},
            'priority_boundary': 'Selected restricted chemistry-only certificate membership; not necessity or guaranteed gain.'})
    return rows, audit


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-purine-precursor-audit', *PRIOR,
             'phase1-new-references', 'phase1-route-references', 'phase1-catalog-references')]
    reports = {str(p): json.loads(p.read_text()) for p in paths}
    directions = Path('data/raw/phase1-reference-discovery/rhea-directions.tsv')
    paths.append(directions)
    families = direction_families(directions.read_text())
    source = reports[str(paths[0])]
    prior = {str(Path('data/reports', n + '.json')): reports[str(Path('data/reports', n + '.json'))] for n in PRIOR}
    rows, audit = queue(source, prior, families)
    masters = {f['RHEA_ID_MASTER'] for r in rows for f in r['rhea_families'].values()}
    lookups, covered, urls = [], set(), set()
    for name in ('phase1-new-references', 'phase1-route-references', 'phase1-catalog-references'):
        filename = str(Path('data/reports', name + '.json'))
        for lookup in reports[filename]['lookups']:
            if lookup['status'] != 'retrieved' or not masters.intersection(lookup['requested_master_ids']) or lookup['url'] in urls:
                continue
            lookups.append({**lookup, 'reused_from': filename})
            covered.update(lookup['requested_master_ids']); urls.add(lookup['url'])
    missing = sorted(masters - covered)
    raw = Path('data/raw/phase1-purine-gap-references'); raw.mkdir(parents=True, exist_ok=True)
    for start in range(0, len(missing), 25):
        batch = missing[start:start + 25]
        expression = '(' + ' OR '.join(f'cc_catalytic_activity:"{rid.lower()}"' for rid in batch) + ') AND reviewed:true AND fragment:false'
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query': expression, 'format': 'tsv', 'fields': 'accession,rhea,organism_name,protein_name'})
        cache = raw / (hashlib.sha256(url.encode()).hexdigest() + '.json')
        if cache.exists():
            lookup = json.loads(cache.read_text())
            if lookup['url'] != url or lookup['requested_master_ids'] != batch:
                raise ValueError('Request cache mismatch')
        else:
            lookup = {'url': url, 'requested_master_ids': batch, 'retrieved_at': datetime.now(timezone.utc).isoformat()}
            try:
                with urllib.request.urlopen(url, timeout=45) as response:
                    payload = response.read()
                exact_annotations(payload.decode(), set(batch))
                sha = hashlib.sha256(payload).hexdigest()
                snapshot = raw / (sha + '.tsv'); snapshot.write_bytes(payload)
                lookup.update(status='retrieved', snapshot=str(snapshot), sha256=sha)
            except (OSError, ValueError) as error:
                lookup.update(status='retrieval-failed', reason=str(error))
            cache.write_text(json.dumps(lookup, separators=(',', ':')) + '\n')
        lookups.append(lookup)
    proteins = attach(rows, lookups)
    for row in rows:
        if not row['rhea_families']:
            row['lookup_status'] = 'no-published-Rhea-family-mapping'
    used = {m['accession'] for r in rows for m in r['reference_matches']}
    report = {'schema': 'cannabis-carbon.phase1-purine-gap-references.v1', 'rows': rows, 'prior_screen_audit': audit,
        'reference_proteins': [proteins[a] for a in sorted(used)], 'lookups': lookups,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'candidate_gap_equations': len(audit), 'previously_screened_equations': len(audit) - len(rows),
            'new_equation_gaps': len(rows), 'requested_master_families': len(masters),
            'reused_lookups': sum('reused_from' in l for l in lookups), 'new_lookup_batches': sum('reused_from' not in l for l in lookups),
            'failed_lookups': sum(l['status'] != 'retrieved' for l in lookups),
            'equations_with_reference_leads': sum(bool(r['reference_matches']) for r in rows),
            'distinct_reference_proteins': len(used), 'lookup_status_counts': dict(Counter(r['lookup_status'] for r in rows))},
        'claim_boundary': 'Reviewed reference annotations, not characterized Cannabis enzymes. All 70 selected gaps retain prior search evidence. Missing reference sequences, weak hits and no hits are distinct from biological absence. No direction, chemistry, precursor supply or completeness claims are promoted. Atom tracing remains deferred.'}
    payload = json.dumps(report, separators=(',', ':')) + '\n'; sha = hashlib.sha256(payload.encode()).hexdigest()
    Path('data/reports/phase1-purine-gap-references.json').write_text(payload)
    groups = [('equation_gap', 'rows', 'reaction_id'), ('prior_screen_audit', 'prior_screen_audit', 'reaction_id'), ('reference_protein', 'reference_proteins', 'accession'), ('lookup', 'lookups', 'url')]
    records = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, key, id_key in groups:
        records.extend((kind, r[id_key], r) for r in report[key])
    count = write_rows(records, sha, Path('data/derived/phase1-purine-gap-references.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
