"""Distinguish unqueried reference families from completed negative searches."""
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


def queue(catalog, supplement, screens, lookups, families):
    added = {e['reaction_id'] for e in supplement['enzyme_evidence']}
    reactions = {r['id']: r for r in catalog['reactions']}
    indexes = {name: {r['reaction_id']: r for r in report['rows']} for name, report in screens.items()}
    covered = {m for l in lookups if l['status'] == 'retrieved' for m in l['requested_master_ids']}
    rows, audit = [], []
    for gap in catalog['gap_priorities']:
        if gap['id'] in added:
            continue
        prior = [{'report': name, **index[gap['id']]} for name, index in indexes.items() if gap['id'] in index]
        if not prior or any(r['passing_alignment_ids'] for r in prior):
            raise ValueError('Remaining gap/prior screen mismatch')
        sources = reactions[gap['id']]['sources']
        ids = sorted({s['source_reaction_id'] for s in sources})
        mapped = {sid: families[sid] for sid in ids if sid in families}
        masters = {f['RHEA_ID_MASTER'] for f in mapped.values()}
        present = {acc for p in prior for acc in p['reference_sequences_present']}
        status = 'reference-sequences-already-screened' if present else 'reference-family-not-queried' if masters-covered else 'reviewed-family-lookup-completed' if masters else 'no-published-family-mapping'
        audit.append({**gap, 'reference_gap_status': status,
            'prior_screens': [{k:p[k] for k in ('report', 'search_status', 'reference_sequences_present', 'reference_sequences_missing', 'raw_alignment_count')} for p in prior],
            'rhea_families': mapped, 'unqueried_master_ids': sorted(masters-covered)})
        if present:
            continue
        r = reactions[gap['id']]
        rows.append({'reaction_id': r['id'], 'left': r['left'], 'right': r['right'],
            'sources': sources, 'source_reaction_ids': ids, 'rhea_families': mapped,
            'target_ids': gap['selected_certificate_target_ids'], 'priority_target_ids': gap['selected_certificate_target_ids'],
            'hypothesis_ids': [], 'previous_reference_gap_status': status,
            'claim_boundary': 'Backfill missing reference discovery, not a repeat of a completed sequence alignment screen.'})
    return rows, audit


def run():
    names = ['phase1-catalog-net-gaps', 'phase1-catalog-evidence',
        'phase1-new-protein-search', 'phase1-route-protein-search', 'phase1-catalog-protein-search',
        'phase1-new-references', 'phase1-route-references', 'phase1-catalog-references']
    paths = [Path('data/reports', n+'.json') for n in names]
    reports = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for report in reports:
        for path, digest in report.get('source_sha256', {}).items():
            if path in hashes and hashes[path] != digest:
                raise ValueError('Backfill source lineage mismatch')
        if report.get('source_discovery') in hashes and report['source_discovery_sha256'] != hashes[report['source_discovery']]:
            raise ValueError('Backfill search/discovery lineage mismatch')
    direction_path = Path(reports[-1]['rhea_direction_source']['snapshot'])
    raw_direction = direction_path.read_bytes()
    if hashlib.sha256(raw_direction).hexdigest() != reports[-1]['rhea_direction_source']['sha256']:
        raise ValueError('Direction snapshot mismatch')
    all_lookups, seen = [], set()
    for i in (5, 6, 7):
        for lookup in reports[i]['lookups']:
            if lookup['url'] in seen:
                continue
            seen.add(lookup['url']); all_lookups.append({**lookup, 'reused_from': str(paths[i])})
    rows, audit = queue(reports[0], reports[1], {str(paths[i]):reports[i] for i in (2,3,4)}, all_lookups, direction_families(raw_direction.decode()))
    masters = {f['RHEA_ID_MASTER'] for r in rows for f in r['rhea_families'].values()}
    lookups = [l for l in all_lookups if masters.intersection(l['requested_master_ids']) and l['status']=='retrieved']
    for l in lookups:
        if hashlib.sha256(Path(l['snapshot']).read_bytes()).hexdigest() != l['sha256']:
            raise ValueError('Cached reference lookup checksum mismatch')
    covered = {m for l in lookups for m in l['requested_master_ids']}
    missing = sorted(masters-covered)
    raw = Path('data/raw/phase1-reference-backfill'); raw.mkdir(parents=True, exist_ok=True)

    def lookup(batch):
        expression = '(' + ' OR '.join(f'cc_catalytic_activity:"{m.lower()}"' for m in batch) + ') AND reviewed:true AND fragment:false'
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query':expression, 'format':'tsv', 'fields':'accession,rhea,organism_name,protein_name'})
        cache = raw / (hashlib.sha256(url.encode()).hexdigest()+'.json')
        if cache.exists():
            result = json.loads(cache.read_text())
            if result['url'] != url or result['requested_master_ids'] != batch:
                raise ValueError('Backfill cache request mismatch')
            return result
        result = {'requested_master_ids':batch, 'url':url, 'retrieved_at':datetime.now(timezone.utc).isoformat()}
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                data = response.read(); result['http_status'] = response.status
            exact_annotations(data.decode(), set(batch))
            digest = hashlib.sha256(data).hexdigest(); snapshot = raw / (digest+'.tsv'); snapshot.write_bytes(data)
            result.update(status='retrieved', snapshot=str(snapshot), sha256=digest)
        except (OSError, ValueError) as error:
            result.update(status='retrieval-failed', reason=str(error))
        cache.write_text(json.dumps(result,separators=(',',':'))+'\n')
        return result

    batches = [missing[i:i+25] for i in range(0,len(missing),25)]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for i, result in enumerate(pool.map(lookup,batches),1):
            lookups.append(result); print(f'Backfill lookup {i}/{len(batches)}: {result["status"]}',flush=True)
    proteins = attach(rows,lookups)
    for r in rows:
        if not r['rhea_families']:
            r['lookup_status'] = 'no-published-family-mapping'
    used = {m['accession'] for r in rows for m in r['reference_matches']}
    result = {'schema':'cannabis-carbon.phase1-reference-backfill.v1', 'rows':rows,
        'remaining_gap_audit':audit, 'reference_proteins':[proteins[a] for a in sorted(used)], 'lookups':lookups,
        'rhea_direction_source':reports[-1]['rhea_direction_source'],
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths+[direction_path]},
        'summary':{'remaining_gap_equations':len(audit), 'reference_gap_status_counts':dict(Counter(r['reference_gap_status'] for r in audit)),
            'no_sequence_equations':len(rows), 'new_master_queries':len(missing), 'new_lookup_batches':len(batches),
            'failed_lookups':sum(l['status']!='retrieved' for l in lookups),
            'equations_with_reference_leads':sum(bool(r['reference_matches']) for r in rows),
            'reference_proteins':len(used), 'lookup_status_counts':dict(Counter(r['lookup_status'] for r in rows))},
        'claim_boundary':'No-reference-sequence is not necessarily a completed reference lookup. This backfill queries previously unqueried published reaction families and reuses checksummed completed lookups. Reviewed annotations are reference leads, not demonstrated activity or physiological direction. Prior failed/weak/negative alignments and frozen chemistry remain unchanged. Atom tracing remains deferred.'}
    payload = json.dumps(result,separators=(',',':'))+'\n'
    Path('data/reports/phase1-reference-backfill.json').write_text(payload)
    groups = [('equation_gap','rows','reaction_id'),('gap_audit','remaining_gap_audit','id'),('reference_protein','reference_proteins','accession'),('lookup','lookups','url')]
    records = [('metadata','report',{k:v for k,v in result.items() if k not in {g[1] for g in groups}})]
    for kind, collection, key in groups:
        records.extend((kind,r[key],r) for r in result[collection])
    count = write_rows(records,hashlib.sha256(payload.encode()).hexdigest(),Path('data/derived/phase1-reference-backfill.ndjson'))
    print(json.dumps({'summary':result['summary'],'rows':count}),flush=True)


if __name__ == '__main__':
    run()
