"""Exact-family reference discovery for previously unscreened catalog net gaps."""
import hashlib
import json
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from .phase1_new_references import attach
from .phase1_reference_discovery import direction_families, exact_annotations
from .phase1_scope import write_rows


def queue(catalog, prior, families):
    reactions = {r['id']: r for r in catalog['reactions']}
    indexes = {name: {r['reaction_id']: r for r in report['rows']} for name, report in prior.items()}
    audit, rows = [], []
    for gap in catalog['gap_priorities']:
        rid = gap['reaction_id']
        reaction = reactions[rid]
        if reaction['enzyme_evidence_ids'] or not reaction['missing_candidate_evidence']:
            raise ValueError('Gap already has candidate evidence')
        previous = [{'report': name, 'reaction_id': rid, 'search_status': index[rid]['search_status'],
            'reference_sequences_present': index[rid]['reference_sequences_present'],
            'passing_alignment_ids': index[rid]['passing_alignment_ids']} for name, index in indexes.items() if rid in index]
        if any(p['passing_alignment_ids'] for p in previous):
            raise ValueError('Prior passing candidate not accounted for in catalog baseline')
        audit.append({**gap, 'prior_screens': previous,
            'disposition': 'retain-prior-search-result' if previous else 'new-reference-discovery'})
        if previous:
            continue
        source_ids = sorted({s['source_reaction_id'].upper() for s in reaction['sources']})
        rows.append({'reaction_id': rid, 'left': reaction['left'], 'right': reaction['right'],
            'sources': reaction['sources'], 'source_reaction_ids': source_ids,
            'target_ids': gap['selected_certificate_target_ids'],
            'priority_target_ids': gap['selected_certificate_target_ids'], 'hypothesis_ids': [],
            'selected_certificate_target_count': gap['selected_certificate_target_count'],
            'rhea_families': {sid: families[sid] for sid in source_ids if sid in families},
            'priority_boundary': 'Selected chemistry-only net certificate membership, not necessity or a guaranteed pathway gain.'})
    return rows, audit


def run():
    names = ['phase1-catalog-net-gaps', 'phase1-route-protein-search', 'phase1-new-protein-search',
        'phase1-reference-discovery', 'phase1-route-references', 'phase1-new-references']
    paths = [Path('data/reports', n + '.json') for n in names]
    reports = [json.loads(p.read_text()) for p in paths]
    for report in reports:
        for filename, digest in report.get('source_sha256', {}).items():
            if hashlib.sha256(Path(filename).read_bytes()).hexdigest() != digest:
                raise ValueError('Pinned source mismatch: ' + filename)
    metadata = reports[3]['rhea_direction_source']
    directions_path = Path(metadata['snapshot'])
    directions = directions_path.read_bytes()
    if hashlib.sha256(directions).hexdigest() != metadata['sha256']:
        raise ValueError('Published direction-family snapshot mismatch')
    rows, audit = queue(reports[0], {str(paths[i]): reports[i] for i in (1, 2)}, direction_families(directions.decode()))
    masters = {f['RHEA_ID_MASTER'] for r in rows for f in r['rhea_families'].values()}
    lookups, covered, urls = [], set(), set()
    for i in (4, 5):
        for lookup in reports[i]['lookups']:
            if lookup['status'] != 'retrieved' or not masters.intersection(lookup['requested_master_ids']) or lookup['url'] in urls:
                continue
            if hashlib.sha256(Path(lookup['snapshot']).read_bytes()).hexdigest() != lookup['sha256']:
                raise ValueError('Cached lookup mismatch')
            lookups.append({**lookup, 'reused_from': str(paths[i])})
            covered.update(lookup['requested_master_ids']); urls.add(lookup['url'])
    missing = sorted(masters - covered)
    raw = Path('data/raw/phase1-catalog-references'); raw.mkdir(parents=True, exist_ok=True)

    def lookup(batch):
        expression = '(' + ' OR '.join(f'cc_catalytic_activity:"{rid.lower()}"' for rid in batch) + ') AND reviewed:true AND fragment:false'
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({
            'query': expression, 'format': 'tsv', 'fields': 'accession,rhea,organism_name,protein_name'})
        cache = raw / (hashlib.sha256(url.encode()).hexdigest() + '.json')
        if cache.exists():
            result = json.loads(cache.read_text())
            if result['url'] != url or result['requested_master_ids'] != batch:
                raise ValueError('Request cache identity mismatch')
            if result['status'] == 'retrieved' and hashlib.sha256(Path(result['snapshot']).read_bytes()).hexdigest() != result['sha256']:
                raise ValueError('Request cache content mismatch')
            return result  # Failed retrievals also remain explicit; no silent retry.
        result = {'requested_master_ids': batch, 'url': url, 'retrieved_at': datetime.now(timezone.utc).isoformat()}
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                data = response.read()
                result['http_status'] = response.status
            exact_annotations(data.decode(), set(batch))
            digest = hashlib.sha256(data).hexdigest()
            snapshot = raw / (digest + '.tsv'); snapshot.write_bytes(data)
            result.update(status='retrieved', sha256=digest, snapshot=str(snapshot))
        except (OSError, ValueError) as error:
            result.update(status='retrieval-failed', reason=str(error))
        cache.write_text(json.dumps(result, separators=(',', ':')) + '\n')
        return result

    batches = [missing[i:i+25] for i in range(0, len(missing), 25)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for i, result in enumerate(pool.map(lookup, batches), 1):
            lookups.append(result)
            print(f'Catalog lookup {i}/{len(batches)}: {result["status"]}', flush=True)
    proteins = attach(rows, lookups)
    for row in rows:
        if not row['rhea_families']:
            row['lookup_status'] = 'no-published-Rhea-family-mapping'
    used = {m['accession'] for row in rows for m in row['reference_matches']}
    output = {'schema': 'cannabis-carbon.phase1-catalog-references.v1', 'rows': rows,
        'prior_screen_audit': audit, 'reference_proteins': [proteins[a] for a in sorted(used)],
        'lookups': lookups, 'rhea_direction_source': metadata,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths + [directions_path]},
        'summary': {'catalog_gap_equations': len(audit), 'previously_screened_equations': len(audit)-len(rows),
            'new_equation_gaps': len(rows), 'requested_master_families': len(masters),
            'reused_lookups': sum('reused_from' in l for l in lookups), 'new_lookup_batches': len(batches),
            'failed_lookups': sum(l['status'] != 'retrieved' for l in lookups),
            'equations_with_reference_leads': sum(bool(r['reference_matches']) for r in rows),
            'distinct_reference_proteins': len(used), 'lookup_status_counts': dict(Counter(r['lookup_status'] for r in rows))},
        'claim_boundary': 'Reviewed exact published reaction-family annotations are reference leads, not Cannabis activity, exact specificity or physiological direction. Prior no-hit, weak-hit and no-reference results remain explicit. Selected-certificate membership is not necessity. Atom tracing remains deferred.'}
    payload = json.dumps(output, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-catalog-references.json').write_text(payload)
    groups = [('equation_gap', 'rows', 'reaction_id'), ('prior_screen_audit', 'prior_screen_audit', 'id'),
        ('reference_protein', 'reference_proteins', 'accession'), ('lookup', 'lookups', 'url')]
    records = [('metadata', 'report', {k:v for k,v in output.items() if k not in {g[1] for g in groups}})]
    for kind, collection, key in groups:
        records.extend((kind, row[key], row) for row in output[collection])
    count = write_rows(records, hashlib.sha256(payload.encode()).hexdigest(), Path('data/derived/phase1-catalog-references.ndjson'))
    print(json.dumps({'summary': output['summary'], 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
