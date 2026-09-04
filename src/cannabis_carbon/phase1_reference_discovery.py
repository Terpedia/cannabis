"""Discover exact Rhea reference annotations for Phase 1 enzyme gaps."""
import csv
import hashlib
import io
import json
import os
import re
import subprocess
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def exact_annotations(tsv, requested):
    """The query can be broad; only explicit returned Rhea IDs are joined."""
    matches = {rid: [] for rid in requested}
    reader = csv.DictReader(io.StringIO(tsv), delimiter='\t')
    if not {'Entry', 'Rhea ID', 'Organism', 'Protein names'} <= set(reader.fieldnames or []):
        raise ValueError('Unexpected UniProt TSV schema')
    for row in reader:
        ids = set(re.findall(r'RHEA:\d+', row['Rhea ID'], re.I))
        ids = {rid.upper() for rid in ids}
        for rid in ids & matches.keys():
            matches[rid].append({
                'accession': row['Entry'], 'annotated_rhea_ids': sorted(ids),
                'organism': row['Organism'], 'protein_name': row['Protein names'],
                'source_url': f"https://www.uniprot.org/uniprotkb/{row['Entry']}/entry",
                'evidence_status': 'reviewed-reference-annotation-not-Cannabis-activity',
            })
    return matches


def direction_families(tsv):
    reader = csv.DictReader(io.StringIO(tsv), delimiter='\t')
    columns = ['RHEA_ID_MASTER', 'RHEA_ID_LR', 'RHEA_ID_RL', 'RHEA_ID_BI']
    if not set(columns) <= set(reader.fieldnames or []):
        raise ValueError('Unexpected Rhea directions schema')
    families = {}
    for row in reader:
        family = {column: 'RHEA:' + row[column] for column in columns}
        for rid in family.values():
            if rid in families:
                raise ValueError('Ambiguous Rhea direction family')
            families[rid] = family
    return families


def attach_families(report, directions):
    families = direction_families(directions)
    proteins = {}
    for lookup in report['uniprot_lookups']:
        if lookup['status'] != 'retrieved':
            continue
        data = Path(lookup['snapshot']).read_bytes()
        if hashlib.sha256(data).hexdigest() != lookup['sha256']:
            raise ValueError('UniProt snapshot checksum mismatch')
        returned_ids = set(re.findall(r'RHEA:\d+', data.decode(), re.I))
        for matches in exact_annotations(data.decode(), returned_ids).values():
            for match in matches:
                record = proteins.setdefault(match['accession'], {**match, 'retrieval_evidence': []})
                provenance = {'url': lookup['url'], 'sha256': lookup['sha256'], 'snapshot': lookup['snapshot']}
                if provenance not in record['retrieval_evidence']:
                    record['retrieval_evidence'].append(provenance)
    for row in report['rows']:
        family = families.get(row['reaction_id'], {})
        row['rhea_direction_family'] = family
        row['family_reference_annotations'] = [
            {**protein, 'join_method': 'explicit-published-Rhea-direction-family',
             'matched_family_rhea_ids': sorted(set(protein['annotated_rhea_ids']) & set(family.values())),
             'direction_status': 'not-established-for-requested-direction'}
            for protein in proteins.values()
            if set(protein['annotated_rhea_ids']) & set(family.values())
            and row['reaction_id'] not in protein['annotated_rhea_ids']]
    report['summary']['variants_with_direction_family_reference_candidates'] = sum(bool(r['family_reference_annotations']) for r in report['rows'])
    report['summary']['distinct_direction_family_reference_proteins'] = len({p['accession'] for r in report['rows'] for p in r['family_reference_annotations']})
    return report


def enrich_saved_report():
    path = Path('data/reports/phase1-reference-discovery.json')
    report = json.loads(path.read_text())
    url = 'https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv'
    with urllib.request.urlopen(url, timeout=45) as response:
        data = response.read()
    snapshot = Path('data/raw/phase1-reference-discovery/rhea-directions.tsv')
    snapshot.write_bytes(data)
    report = attach_families(report, data.decode())
    report['rhea_direction_source'] = {'url': url, 'snapshot': str(snapshot), 'sha256': hashlib.sha256(data).hexdigest()}
    path.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(report['summary'])


def export_table(report_path, output):
    """One row per reaction variant, retaining reference provenance as JSON."""
    data = report_path.read_bytes()
    report = json.loads(data)
    digest = hashlib.sha256(data).hexdigest()
    rows = []
    for row in report['rows']:
        rows.append({
            'reaction_id': row['reaction_id'], 'reaction_smarts': row['reaction_smarts'],
            'exact_reference_status': row['status'],
            'family_reference_count': len(row['family_reference_annotations']),
            'sequence_search_status': row['sequence_search_status'],
            'exact_references_json': json.dumps(row['reference_annotations']),
            'family_references_json': json.dumps(row['family_reference_annotations']),
            'direction_family_json': json.dumps(row['rhea_direction_family']),
            'terpedia_catalog_rows_json': json.dumps(row['terpedia_exact_variant_rows']),
            'report_sha256': digest, 'source_overlay_sha256': report['source_overlay_sha256'],
            'direction_mapping_sha256': report['rhea_direction_source']['sha256'],
            'claim_boundary': row['claim_boundary'],
        })
    output.write_text(''.join(json.dumps(row, separators=(',', ':')) + '\n' for row in rows))
    return len(rows)


def run():
    source = Path('data/reports/phase1-map-evidence.json')
    gaps = [r for r in json.loads(source.read_text())['rows'] if r['evidence_status'] == 'missing-reference']
    ids = sorted({r['reaction_id'] for r in gaps})
    if not ids or not all(re.fullmatch(r'RHEA:\d+', rid) for rid in ids):
        raise ValueError('Expected explicit Rhea gap identifiers')
    raw = Path('data/raw/phase1-reference-discovery')
    raw.mkdir(parents=True, exist_ok=True)
    table = 'terpedia-489015.terpedia_core.terpene_reaction_smarts_catalog_normalized_current_v2'
    query = f"SELECT rule_id, reaction_smarts, original_reaction_smarts, source_uniprot_id, source_genbank_id, source_ec_number, direction_mode, source_url, source_evidence_type FROM `{table}` WHERE rule_id IN (" + ','.join("'" + rid + "'" for rid in ids) + ')'
    catalog_bytes = subprocess.check_output([
        os.environ.get('CANNABIS_BQ', 'bq'), '--format=json', 'query', '--use_legacy_sql=false',
        '--maximum_bytes_billed=1073741824', '--max_rows=10000', query])
    catalog = json.loads(catalog_bytes)
    catalog_path = raw / 'terpedia-catalog.json'
    catalog_path.write_bytes(catalog_bytes)
    lookups, index = [], {}
    for start in range(0, len(ids), 25):
        batch = ids[start:start + 25]
        expression = '(' + ' OR '.join(f'cc_catalytic_activity:"{rid.lower()}"' for rid in batch) + ') AND reviewed:true AND fragment:false'
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({
            'query': expression, 'format': 'tsv', 'fields': 'accession,rhea,organism_name,protein_name'})
        lookup = {'requested_reaction_ids': batch, 'url': url}
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                data = response.read()
            parsed = exact_annotations(data.decode(), set(batch))
            digest = hashlib.sha256(data).hexdigest()
            snapshot = raw / (digest + '.tsv')
            snapshot.write_bytes(data)
            lookup.update(status='retrieved', sha256=digest, snapshot=str(snapshot))
            index.update(parsed)
        except (OSError, ValueError) as error:
            lookup.update(status='retrieval-failed', reason=str(error))
        lookups.append(lookup)
        print(f'Reference lookup {min(start + 25, len(ids))}/{len(ids)}: {lookup["status"]}', flush=True)
    rows = []
    for gap in gaps:
        rid = gap['reaction_id']
        matches = index.get(rid, [])
        exact_catalog = [r for r in catalog if r['rule_id'] == rid and r['reaction_smarts'] == gap['reaction_smarts']]
        rows.append({'reaction_id': rid, 'reaction_smarts': gap['reaction_smarts'],
                     'terpedia_exact_variant_rows': exact_catalog,
                     'reference_annotations': matches,
                     'status': 'reference-annotations-found' if matches else 'no-reviewed-exact-reference' if rid in index else 'retrieval-failed',
                     'sequence_search_status': 'not-yet-searched',
                     'claim_boundary': 'Exact Rhea annotation identifies a reference-search input, not proof of stereo-insensitive endpoint specificity, physiological direction, or Cannabis activity.'})
    result = {'schema': 'cannabis-carbon.phase1-reference-discovery.v1',
              'generated_at': datetime.now(timezone.utc).isoformat(),
              'source_overlay': str(source), 'source_overlay_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'terpedia_catalog': {'table': table, 'query': query, 'snapshot': str(catalog_path),
                                  'sha256': hashlib.sha256(catalog_bytes).hexdigest(), 'row_count': len(catalog)},
              'uniprot_lookups': lookups,
              'summary': {'gap_variants': len(rows),
                          'variants_with_exact_terpedia_catalog_match': sum(bool(r['terpedia_exact_variant_rows']) for r in rows),
                          'variants_with_reviewed_reference_annotations': sum(bool(r['reference_annotations']) for r in rows),
                          'distinct_reference_proteins': len({p['accession'] for r in rows for p in r['reference_annotations']}),
                          'failed_lookup_variants': sum(r['status'] == 'retrieval-failed' for r in rows)},
              'rows': rows}
    Path('data/reports/phase1-reference-discovery.json').write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(result['summary'])


if __name__ == '__main__':
    run()
    enrich_saved_report()
