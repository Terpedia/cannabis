"""Provenance-only Rhea annotation triage of remaining catalog enzyme gaps.

EC links establish a source classification, not Cannabis activity. Missing EC
links never establish spontaneity. Full source triples remain available for
curation; this report does not change equations, evidence or certificates.
"""
import hashlib
import json
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .phase1_reference_discovery import direction_families
from .phase1_scope import write_rows

RH = 'http://rdf.rhea-db.org/'
ENDPOINT = 'https://sparql.rhea-db.org/sparql'
CATALOG = Path('data/reports/phase1-catalog-net-gaps.json')
EVIDENCE = Path('data/reports/phase1-combined-catalog-evidence.json')
DIRECTIONS = Path('data/raw/phase1-reference-discovery/rhea-directions.tsv')
RAW = Path('data/raw/phase1-gap-annotations')


def queue(catalog, evidence, families):
    reactions = {r['id']: r for r in catalog['reactions']}
    added = {e['reaction_id'] for e in evidence['enzyme_evidence']}
    rows = []
    for gap in catalog['gap_priorities']:
        rid = gap['reaction_id']
        if rid in added:
            continue
        reaction = reactions[rid]
        joins = []
        for source in reaction['sources']:
            sid = source['source_reaction_id']
            family = families.get(sid)
            joins.append({'source': source, 'published_direction_family': family,
                          'master_id': family['RHEA_ID_MASTER'] if family else None})
        rows.append({**gap, 'source_joins': joins})
    if len(rows) != evidence['summary']['remaining_missing_candidate_equations']:
        raise ValueError('Remaining gap denominator mismatch')
    return rows


def fetch(masters):
    """Save each successful response and request metadata; cached bytes are pinned."""
    RAW.mkdir(parents=True, exist_ok=True)
    lookups = []
    for start in range(0, len(masters), 25):
        batch = masters[start:start + 25]
        values = ' '.join('<' + RH + m.split(':')[1] + '>' for m in batch)
        query = f'SELECT ?s ?p ?o WHERE {{ VALUES ?s {{ {values} }} ?s ?p ?o }} ORDER BY ?s ?p ?o'
        key = hashlib.sha256(query.encode()).hexdigest()
        metadata = RAW / (key + '.request.json')
        snapshot = RAW / (key + '.response.json')
        if metadata.exists():
            lookup = json.loads(metadata.read_text())
            if lookup['query'] != query or lookup['requested_master_ids'] != batch:
                raise ValueError('Cached request mismatch')
            if hashlib.sha256(Path(lookup['snapshot']).read_bytes()).hexdigest() != lookup['sha256']:
                raise ValueError('Rhea snapshot checksum mismatch')
        else:
            url = ENDPOINT + '?' + urllib.parse.urlencode({'query': query, 'format': 'json'})
            with urllib.request.urlopen(url, timeout=45) as response:
                payload = response.read()
            data = json.loads(payload)
            if data['head']['vars'] != ['s', 'p', 'o']:
                raise ValueError('Unexpected SPARQL result schema')
            snapshot.write_bytes(payload)
            lookup = {'query': query, 'url': url, 'requested_master_ids': batch,
                      'snapshot': str(snapshot), 'sha256': hashlib.sha256(payload).hexdigest(),
                      'retrieved_at': datetime.now(timezone.utc).isoformat(), 'status': 'retrieved'}
            metadata.write_text(json.dumps(lookup, indent=2) + '\n')
        lookups.append(lookup)
    return lookups


def assemble(rows, lookups):
    annotations = {}
    for lookup in lookups:
        payload = Path(lookup['snapshot']).read_bytes()
        if hashlib.sha256(payload).hexdigest() != lookup['sha256']:
            raise ValueError('Rhea snapshot checksum mismatch')
        requested = set(lookup['requested_master_ids'])
        for binding in json.loads(payload)['results']['bindings']:
            subject = binding['s']
            if subject['type'] != 'uri' or not subject['value'].startswith(RH):
                raise ValueError('Unexpected Rhea subject')
            master = 'RHEA:' + subject['value'][len(RH):]
            if master not in requested:
                raise ValueError('Unrequested Rhea subject')
            record = annotations.setdefault(master, {'id': master, 'triples': [],
                'source_url': 'https://www.rhea-db.org/rhea/' + master.split(':')[1],
                'retrieval': lookup})
            record['triples'].append({'predicate': binding['p'], 'object': binding['o']})

    def objects(master, predicate):
        return sorted({t['object']['value'] for t in annotations.get(master, {}).get('triples', [])
                       if t['predicate']['value'] == predicate})

    output = []
    for row in rows:
        masters = sorted({j['master_id'] for j in row['source_joins'] if j['master_id']})
        ec = sorted({v for m in masters for v in objects(m, RH + 'ec')})
        citations = sorted({v for m in masters for v in objects(m, RH + 'citation')})
        missing = [m for m in masters if m not in annotations]
        unmapped = sorted({j['source']['source_reaction_id'] for j in row['source_joins'] if not j['master_id']})
        status = ('source-ec-linked-Cannabis-function-unresolved' if ec else
                  'source-annotation-incomplete' if missing or unmapped else
                  'no-source-ec-link-catalysis-unresolved')
        output.append({**row, 'annotation_master_ids': masters, 'ec_annotation_urls': ec,
            'citation_urls': citations, 'source_equations': {m: objects(m, RH + 'equation') for m in masters},
            'missing_annotation_master_ids': missing, 'unmapped_source_reaction_ids': unmapped,
            'annotation_status': status, 'spontaneous_status': 'unresolved-not-inferred-from-missing-EC-or-proteins',
            'next_curation_test': 'Read linked primary literature for organism, enzyme complex, cofactors and measured direction. Check plant-specific alternatives before prioritizing genome searches. No EC link is not evidence of spontaneous chemistry.'})
    return {'schema': 'cannabis-carbon.phase1-gap-annotations.v1', 'rows': output,
            'source_annotations': [annotations[k] for k in sorted(annotations)], 'lookups': lookups,
            'summary': {'remaining_gap_equations': len(output),
                        'source_master_annotations': len(annotations),
                        'equations_with_citations': sum(bool(r['citation_urls']) for r in output),
                        'annotation_status_counts': dict(Counter(r['annotation_status'] for r in output)),
                        'new_Cannabis_activity_claims': 0, 'new_spontaneous_claims': 0},
            'claim_boundary': 'Source annotation triage only. EC links and publications do not establish Cannabis activity, physiological direction, spontaneous chemistry or pathway necessity. Ranking counts selected certificates, not essential reactions. Equations, candidate evidence and net certificates are unchanged.'}


def run():
    catalog, evidence = [json.loads(p.read_text()) for p in (CATALOG, EVIDENCE)]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (CATALOG, EVIDENCE, DIRECTIONS)}
    if evidence['source_sha256'][str(CATALOG)] != hashes[str(CATALOG)]:
        raise ValueError('Catalog evidence lineage mismatch')
    rows = queue(catalog, evidence, direction_families(DIRECTIONS.read_text()))
    masters = sorted({j['master_id'] for r in rows for j in r['source_joins'] if j['master_id']})
    report = assemble(rows, fetch(masters))
    report['source_sha256'] = hashes
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-gap-annotations.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    records = [('metadata', 'audit', {k: v for k, v in report.items() if k not in ('rows', 'source_annotations', 'lookups')})]
    records += [('equation_gap', r['id'], r) for r in report['rows']]
    records += [('source_annotation', r['id'], r) for r in report['source_annotations']]
    records += [('lookup', str(i), r) for i, r in enumerate(report['lookups'])]
    count = write_rows(records, sha, Path('data/derived/phase1-gap-annotations.ndjson'))
    print(json.dumps({'summary': report['summary'], 'rows': count, 'sha256': sha}))


if __name__ == '__main__':
    run()
