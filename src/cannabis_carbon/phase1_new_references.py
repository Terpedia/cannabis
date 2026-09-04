"""Reviewed reaction-family reference discovery for carbon-bearing enzyme gaps."""
import hashlib
import json
import re
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from .phase1_reference_discovery import direction_families, exact_annotations


def queue(report, families):
    compounds = {c['id']: c for c in report['compounds']}
    reactions = {r['id']: r for r in report['reactions']}
    supported = {h['cannabisdb_id'] for h in report['hypotheses'] if h['has_candidate_enzyme_evidence']}
    rows = {}
    for h in report['hypotheses']:
        if not compounds[h['compound_id']]['carbon_count'] or h['has_candidate_enzyme_evidence']:
            continue
        rid = h['reaction_id']
        reaction = reactions[rid]
        source_ids = sorted({s['source_reaction_id'].upper() for s in reaction['sources']})
        row = rows.setdefault(rid, {'reaction_id': rid, 'left': reaction['left'], 'right': reaction['right'],
            'source_reaction_ids': source_ids, 'sources': reaction['sources'], 'hypothesis_ids': [],
            'target_ids': [], 'priority_target_ids': [],
            'rhea_families': {sid: families[sid] for sid in source_ids if sid in families}})
        row['hypothesis_ids'].append(h['id'])
        if h['cannabisdb_id'] not in row['target_ids']:
            row['target_ids'].append(h['cannabisdb_id'])
        if h['cannabisdb_id'] not in supported and h['cannabisdb_id'] not in row['priority_target_ids']:
            row['priority_target_ids'].append(h['cannabisdb_id'])
    return list(rows.values())


def attach(rows, lookups):
    proteins = {}
    queried = set()
    successful = set()
    for lookup in lookups:
        queried.update(lookup['requested_master_ids'])
        if lookup['status'] != 'retrieved':
            continue
        successful.update(lookup['requested_master_ids'])
        raw = Path(lookup['snapshot']).read_bytes()
        if hashlib.sha256(raw).hexdigest() != lookup['sha256']:
            raise ValueError('Reference snapshot checksum mismatch')
        ids = {x.upper() for x in re.findall(r'RHEA:\d+', raw.decode(), re.I)}
        for matches in exact_annotations(raw.decode(), ids).values():
            for match in matches:
                protein = proteins.setdefault(match['accession'], {**match, 'retrieval_evidence': []})
                evidence = {k: lookup[k] for k in ('url', 'snapshot', 'sha256')}
                if evidence not in protein['retrieval_evidence']:
                    protein['retrieval_evidence'].append(evidence)
    index = {}
    for protein in proteins.values():
        for rid in protein['annotated_rhea_ids']:
            index.setdefault(rid, set()).add(protein['accession'])
    for row in rows:
        family_ids = {rid for family in row['rhea_families'].values() for rid in family.values()}
        masters = {family['RHEA_ID_MASTER'] for family in row['rhea_families'].values()}
        accessions = sorted({acc for rid in family_ids for acc in index.get(rid, [])})
        row['reference_matches'] = [{'accession': acc,
            'exact_source_id_matches': sorted(set(proteins[acc]['annotated_rhea_ids']) & set(row['source_reaction_ids'])),
            'family_annotation_matches': sorted(set(proteins[acc]['annotated_rhea_ids']) & family_ids),
            'direction_status': 'not-established-for-hypothesis',
            'join_method': 'explicit-published-Rhea-family; no numeric-ID arithmetic'} for acc in accessions]
        row['lookup_status'] = 'references-found' if accessions else 'no-reviewed-reference-returned' if masters and masters <= successful else 'lookup-incomplete-or-failed' if masters & queried else 'not-searched-in-priority-pass'
        row['next_step'] = 'Retrieve validated reference sequences and screen the full Cannabis proteome; evaluate specificity before assay.' if accessions else 'Resolve missing reference evidence; do not interpret missing annotation as biological absence.'
    return proteins


def run():
    report_path = Path('data/reports/phase1-target-hypotheses.json')
    metadata_path = Path('data/reports/phase1-reference-discovery.json')
    metadata = json.loads(metadata_path.read_text())['rhea_direction_source']
    directions_path = Path(metadata['snapshot'])
    directions = directions_path.read_bytes()
    if hashlib.sha256(directions).hexdigest() != metadata['sha256']:
        raise ValueError('Rhea family snapshot checksum mismatch')
    report = json.loads(report_path.read_text())
    rows = queue(report, direction_families(directions.decode()))
    masters = sorted({family['RHEA_ID_MASTER'] for row in rows if row['priority_target_ids'] for family in row['rhea_families'].values()})
    raw = Path('data/raw/phase1-new-references'); raw.mkdir(parents=True, exist_ok=True)
    def lookup(batch):
        expression = '(' + ' OR '.join(f'cc_catalytic_activity:"{rid.lower()}"' for rid in batch) + ') AND reviewed:true AND fragment:false'
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query': expression, 'format': 'tsv', 'fields': 'accession,rhea,organism_name,protein_name'})
        row = {'requested_master_ids': batch, 'url': url}
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                data = response.read()
            exact_annotations(data.decode(), set(batch))
            digest = hashlib.sha256(data).hexdigest()
            path = raw / (digest + '.tsv'); path.write_bytes(data)
            row.update(status='retrieved', sha256=digest, snapshot=str(path))
        except (OSError, ValueError) as error:
            row.update(status='retrieval-failed', reason=str(error))
        return row
    batches = [masters[start:start + 25] for start in range(0, len(masters), 25)]
    lookups = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for n, result in enumerate(pool.map(lookup, batches), 1):
            lookups.append(result)
            print(f'Lookup {n}/{len(batches)}: {result["status"]}', flush=True)
    proteins = attach(rows, lookups)
    output = {'schema': 'cannabis-carbon.phase1-new-references.v1', 'rows': rows,
        'reference_proteins': list(proteins.values()), 'lookups': lookups, 'rhea_direction_source': metadata,
        'summary': {'carbon_bearing_gap_hypotheses': sum(len(r['hypothesis_ids']) for r in rows),
            'balanced_equation_gaps': len(rows), 'priority_targets': len({tid for r in rows for tid in r['priority_target_ids']}),
            'queried_master_families': len(masters), 'failed_batches': sum(l['status'] != 'retrieved' for l in lookups),
            'equations_with_reference_leads': sum(bool(r['reference_matches']) for r in rows),
            'priority_targets_with_reference_leads': len({tid for r in rows if r['reference_matches'] for tid in r['priority_target_ids']}),
            'distinct_reference_proteins': len({m['accession'] for r in rows for m in r['reference_matches']})},
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [report_path, metadata_path, directions_path]},
        'claim_boundary': 'All carbon-bearing hypothesis enzyme gaps retained. Queries prioritize targets with no candidate-supported alternative. Reviewed protein annotations are reference leads from other organisms, not Cannabis activity, exact substrate specificity, physiological direction, or CO2 reachability. Atom tracing remains deferred.'}
    Path('data/reports/phase1-new-references.json').write_text(json.dumps(output, separators=(',', ':')) + '\n')
    print(json.dumps(output['summary']), flush=True)


def export_table(report_path, output):
    raw = report_path.read_bytes()
    report = json.loads(raw)
    digest = hashlib.sha256(raw).hexdigest()
    records = []
    for kind, collection, key in [('equation_gap', 'rows', 'reaction_id'),
                                   ('reference_protein', 'reference_proteins', 'accession'),
                                   ('lookup', 'lookups', 'url')]:
        for row in report[collection]:
            records.append({'record_kind': kind, 'record_id': row[key],
                'status': row.get('lookup_status') or row.get('status'),
                'record_json': json.dumps(row, separators=(',', ':')), 'report_sha256': digest})
    output.write_text(''.join(json.dumps(r, separators=(',', ':')) + '\n' for r in records))
    return len(records)


if __name__ == '__main__':
    run()
