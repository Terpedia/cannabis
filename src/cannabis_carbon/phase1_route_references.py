"""Discover reviewed reference leads for all selected-route enzyme gaps."""
import hashlib
import json
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from .phase1_reference_discovery import direction_families, exact_annotations
from .phase1_new_references import attach
from .phase1_scope import write_rows


def queue(certificates, network, families):
    reactions = {r['id']: r for r in network['reactions']}
    rows = []
    for gap in certificates['enzyme_gap_queue']:
        reaction = reactions[gap['reaction_id']]
        if reaction['enzyme_evidence_ids'] or gap['sources'] != reaction['sources']:
            raise ValueError('Gap/network evidence mismatch')
        source_ids = sorted({s['source_reaction_id'].upper() for s in reaction['sources']})
        rows.append({'reaction_id': reaction['id'], 'left': reaction['left'], 'right': reaction['right'],
            'source_reaction_ids': source_ids, 'sources': reaction['sources'],
            'target_ids': gap['target_ids'], 'priority_target_ids': gap['target_ids'],
            'route_indices': gap['route_indices'], 'hypothesis_ids': [],
            'selected_route_target_count': gap['selected_route_target_count'],
            'rhea_families': {sid: families[sid] for sid in source_ids if sid in families},
            'priority_boundary': 'Selected deterministic route membership, not an unavoidable biological bottleneck.'})
    return rows


def run():
    cert_path = Path('data/reports/phase1-route-certificates.json')
    network_path = Path('data/reports/phase1-full-balanced-network.json')
    metadata_path = Path('data/reports/phase1-reference-discovery.json')
    prior_path = Path('data/reports/phase1-new-references.json')
    cert = json.loads(cert_path.read_text())
    if hashlib.sha256(network_path.read_bytes()).hexdigest() != cert['source_sha256'][str(network_path)]:
        raise ValueError('Certificate network checksum mismatch')
    metadata = json.loads(metadata_path.read_text())['rhea_direction_source']
    directions_path = Path(metadata['snapshot'])
    directions = directions_path.read_bytes()
    if hashlib.sha256(directions).hexdigest() != metadata['sha256']:
        raise ValueError('Rhea direction source mismatch')
    rows = queue(cert, json.loads(network_path.read_text()), direction_families(directions.decode()))
    masters = {family['RHEA_ID_MASTER'] for row in rows for family in row['rhea_families'].values()}
    lookups, covered = [], set()
    for lookup in json.loads(prior_path.read_text())['lookups']:
        if lookup['status'] != 'retrieved' or not masters.intersection(lookup['requested_master_ids']):
            continue
        snapshot = Path(lookup['snapshot'])
        if not snapshot.exists():
            continue
        if hashlib.sha256(snapshot.read_bytes()).hexdigest() != lookup['sha256']:
            raise ValueError('Cached reference lookup checksum mismatch')
        lookups.append({**lookup, 'reused_from': str(prior_path)})
        covered.update(lookup['requested_master_ids'])
    missing = sorted(masters - covered)
    raw = Path('data/raw/phase1-route-references'); raw.mkdir(parents=True, exist_ok=True)

    def lookup(batch):
        expression = '(' + ' OR '.join(f'cc_catalytic_activity:"{rid.lower()}"' for rid in batch) + ') AND reviewed:true AND fragment:false'
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query': expression, 'format': 'tsv', 'fields': 'accession,rhea,organism_name,protein_name'})
        result = {'requested_master_ids': batch, 'url': url}
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                data = response.read()
            exact_annotations(data.decode(), set(batch))
            digest = hashlib.sha256(data).hexdigest()
            snapshot = raw / (digest + '.tsv'); snapshot.write_bytes(data)
            result.update(status='retrieved', sha256=digest, snapshot=str(snapshot))
        except (OSError, ValueError) as error:
            result.update(status='retrieval-failed', reason=str(error))
        return result

    batches = [missing[i:i + 25] for i in range(0, len(missing), 25)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for i, result in enumerate(pool.map(lookup, batches), 1):
            lookups.append(result)
            print(f'Route reference lookup {i}/{len(batches)}: {result["status"]}', flush=True)
    proteins = attach(rows, lookups)
    for row in rows:
        if not row['rhea_families']:
            row['lookup_status'] = 'no-published-Rhea-family-mapping'
    used = {m['accession'] for r in rows for m in r['reference_matches']}
    result = {'schema': 'cannabis-carbon.phase1-route-references.v1', 'rows': rows,
        'reference_proteins': [proteins[acc] for acc in sorted(used)], 'lookups': lookups,
        'rhea_direction_source': metadata,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [cert_path, network_path, metadata_path, prior_path, directions_path]},
        'summary': {'selected_route_equation_gaps': len(rows), 'requested_master_families': len(masters),
            'reused_lookups': sum('reused_from' in l for l in lookups), 'new_lookups': len(batches),
            'failed_lookups': sum(l['status'] != 'retrieved' for l in lookups),
            'equations_with_reference_leads': sum(bool(r['reference_matches']) for r in rows),
            'distinct_reference_proteins': len(used)},
        'claim_boundary': 'Reviewed reaction-family annotations are reference leads, not necessarily direct assays, Cannabis activity, substrate specificity or physiological direction. All selected-route gaps remain explicit. Atom tracing remains deferred.'}
    payload = json.dumps(result, separators=(',', ':')) + '\n'
    for folder in ('data/reports', 'docs/data'):
        Path(folder, 'phase1-route-references.json').write_text(payload)
    records = [('metadata', 'report', {k: v for k, v in result.items() if k not in ('rows', 'reference_proteins', 'lookups')})]
    records.extend(('equation_gap', r['reaction_id'], r) for r in rows)
    records.extend(('reference_protein', p['accession'], p) for p in result['reference_proteins'])
    records.extend(('lookup', l['url'], l) for l in lookups)
    write_rows(records, hashlib.sha256(payload.encode()).hexdigest(), Path('data/reports/phase1-route-references.ndjson'))
    print(json.dumps(result['summary']), flush=True)


if __name__ == '__main__':
    run()
