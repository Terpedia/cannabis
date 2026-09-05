"""Resolve original MARTS archive references without substituting functional annotations."""
import hashlib
import json
import re
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from .phase1_completion_protein_discovery import BOUNDARY
from .phase1_scope import write_rows


def archive_links(discovery, review):
    reviewed = {r['requested_accession']: r for r in review['references']}
    links = []
    for row in discovery['excluded_source_references']:
        original = row['source_uniprot_id']; archive = None; source = None
        if row['status'] == 'UniParc-requires-explicit-resolution':
            archive = original
            method = 'exact-original-MARTS-UniParc-identifier'
        elif row['status'] == 'known-inactive-reference-requires-resolution':
            source = reviewed[original]
            if hashlib.sha256(source['response_text'].encode()).hexdigest() != source['response_sha256'] or json.loads(source['response_text']) != source['entry']:
                raise ValueError('Inactive entry snapshot checksum mismatch')
            entry = source['entry']
            if entry.get('primaryAccession') != original or entry.get('entryType') != 'Inactive':
                raise ValueError('Inactive entry identity mismatch')
            archive = entry.get('extraAttributes', {}).get('uniParcId')
            method = 'explicit-UniProt-inactive-entry-UniParc-link; not an active-accession replacement'
        if archive:
            if not re.fullmatch(r'UPI[0-9A-F]{10}', archive):
                raise ValueError('Invalid explicit archive identifier')
            links.append({**row, 'archive_accession': archive, 'join_method': method,
                'inactive_entry_provenance': {k: source[k] for k in ('requested_accession', 'url', 'response_sha256', 'retrieved_at')} if source else None})
    return links


def validate_entry(accession, entry):
    if entry.get('uniParcId') != accession:
        raise ValueError('Returned UniParc identity differs from request')
    seq = entry.get('sequence', {})
    value = seq.get('value', '')
    if not value or set(value) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO*') or len(value) != seq.get('length'):
        raise ValueError('Invalid archived sequence or length')
    if hashlib.md5(value.encode()).hexdigest().upper() != seq.get('md5', '').upper():
        raise ValueError('Archived sequence MD5 mismatch')
    return {'accession': accession, 'header': accession + ' original-MARTS-UniParc-sequence',
        'sequence': value, 'sequence_sha256': hashlib.sha256(value.encode()).hexdigest(),
        'source_url': 'https://www.uniprot.org/uniparc/' + accession,
        'identity_boundary': 'Exact archived sequence identity only; no inferred active accession, organism-specific function or enzymatic capability.'}


def retrieve(accession, folder):
    path = folder / (accession + '.json')
    if path.exists():
        result = json.loads(path.read_text())
    else:
        url = f'https://rest.uniprot.org/uniparc/{accession}.json'
        result = {'requested_accession': accession, 'url': url,
                  'retrieved_at': datetime.now(timezone.utc).isoformat()}
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                data = response.read(); result['http_status'] = response.status
            result.update(response_text=data.decode(), response_sha256=hashlib.sha256(data).hexdigest(), status='retrieved')
        except (OSError, UnicodeError) as error:
            result.update(status='retrieval-failed', reason=str(error))
        path.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    if result['requested_accession'] != accession:
        raise ValueError('Cached archive request identity mismatch')
    if result['status'] != 'retrieved':
        return result
    if hashlib.sha256(result['response_text'].encode()).hexdigest() != result['response_sha256']:
        raise ValueError('Archive response checksum mismatch')
    try:
        entry = json.loads(result['response_text'])
        reference = validate_entry(accession, entry)
    except (ValueError, KeyError, TypeError) as error:
        return {**result, 'status': 'invalid-archive-response', 'reason': str(error)}
    return {**result, 'reference_sequence': reference, 'cross_references': entry.get('uniParcCrossReferences', [])}


def build(discovery, review, resolved):
    links = archive_links(discovery, review)
    wanted = {r['archive_accession'] for r in links}
    indexed = {r['requested_accession']: r for r in resolved}
    if indexed.keys() != wanted or len(indexed) != len(resolved):
        raise ValueError('Archive resolution inventory mismatch')
    by_completion = defaultdict(list)
    for link in links:
        record = indexed[link['archive_accession']]
        if link['source_uniprot_id'] != link['archive_accession'] and record['status'] == 'retrieved':
            # Require the archived entry itself to independently retain the old
            # accession, not only the inactive entry's forward pointer.
            if not any(x.get('id') == link['source_uniprot_id'] and x.get('database', '').startswith('UniProt') for x in record['cross_references']):
                link = {**link, 'resolution_status': 'inactive-link-not-reciprocal; not screened'}
            else:
                link = {**link, 'resolution_status': 'explicit-reciprocal-archive-link'}
        else:
            link = {**link, 'resolution_status': record['status']}
        by_completion[link['completion_id']].append(link)
    rows = []
    for original in discovery['rows']:
        selected = by_completion.get(original['hypothesis_ids'][0], [])
        if not selected:
            continue
        references = defaultdict(list)
        for link in selected:
            if link['resolution_status'] != 'inactive-link-not-reciprocal; not screened':
                references[link['archive_accession']].append(link)
        rows.append({**original, 'reference_matches': [{'accession': accession,
            'source_record_ids': sorted({r['source_record_id'] for r in matches}),
            'original_source_accessions': sorted({r['source_uniprot_id'] for r in matches}),
            'archive_links': matches, 'join_method': 'original-MARTS-explicit-archive-identity; no template-enzyme transfer',
            'claim_boundary': BOUNDARY} for accession, matches in sorted(references.items())],
            'archive_resolution_links': selected})
    eligible = {m['accession'] for r in rows for m in r['reference_matches']}
    return {'schema': 'cannabis-carbon.phase1-archived-references.v1', 'rows': rows,
        'retrievals': resolved, 'reference_sequences': [r['reference_sequence'] for r in resolved
            if r['status'] == 'retrieved' and r['requested_accession'] in eligible],
        'summary': {'original_source_associations': len(links), 'requested_archive_identifiers': len(wanted),
            'retrieval_status_counts': dict(Counter(r['status'] for r in resolved)),
            'equation_gaps': len(rows), 'source_target_records': len({t for r in rows for t in r['target_ids']}),
            'resolution_link_status_counts': dict(Counter(l['resolution_status'] for r in rows for l in r['archive_resolution_links']))},
        'archive_documentation': ['https://www.uniprot.org/uniparc', 'https://www.uniprot.org/help/api_retrieve_entries'],
        'claim_boundary': 'Archive recovery identifies original source sequences, not enzymes or Cannabis activity. Original UniParc identifiers and explicit reciprocal inactive-entry links are preserved. No current UniProt accession or annotation is substituted. ' + BOUNDARY}


def run():
    paths = [Path('data/reports', n + '.json') for n in ['phase1-completion-protein-discovery', 'phase1-marts-gap-references']]
    discovery, review = [json.loads(p.read_text()) for p in paths]
    accessions = sorted({r['archive_accession'] for r in archive_links(discovery, review)})
    folder = Path('data/raw/phase1-archived-references'); folder.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        resolved = list(pool.map(lambda accession: retrieve(accession, folder), accessions))
    report = build(discovery, review, resolved)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    path = Path('data/reports/phase1-archived-references.json'); path.write_text(payload)
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in ('rows', 'retrievals', 'reference_sequences')})]
    for kind, collection, key in [('equation', 'rows', 'reaction_id'), ('retrieval', 'retrievals', 'requested_accession'), ('reference_sequence', 'reference_sequences', 'accession')]:
        rows.extend((kind, row[key], row) for row in report[collection])
    digest = hashlib.sha256(payload.encode()).hexdigest()
    count = write_rows(rows, digest, Path('data/derived/phase1-archived-references.ndjson'))
    print(json.dumps({'summary': report['summary'], 'rows': count, 'sha256': digest}), flush=True)


if __name__ == '__main__':
    run()
