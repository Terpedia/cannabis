"""Join exact ChEBI endpoint structures after bounded identifier retrieval."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import RDLogger
from .phase1_protonation_source_join import build, exact_key


def run():
    RDLogger.DisableLog('rdApp.warning')
    audit_path = Path('data/reports/phase1-expanded-protonation-audit.json')
    parent_path = Path('data/reports/phase1-protonation-source-join.json')
    ledger_path = Path('data/raw/chebi-protonation-leads-20260906/retrievals.json')
    mapping_path = Path('data/raw/rhea-chebi-ph73-mapping-20260906.tsv')
    structure_path = Path('data/raw/rhea-chebi-smiles-20260906.tsv')
    parent, ledger = [json.loads(p.read_text()) for p in (parent_path, ledger_path)]
    for doc in (parent, ledger):
        for source, sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(source).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed source snapshot')
    expected = {e['source_chebi_id'] for row in parent['rows'] for e in row['unresolved_source_structure_leads']}
    if len(ledger['rows']) != len(expected) or {r['requested_chebi_id'] for r in ledger['rows']} != expected:
        raise ValueError('Incomplete or duplicate retrieval inventory')
    with structure_path.open() as handle:
        structures = list(csv.reader(handle, delimiter='\t'))
    with mapping_path.open() as handle:
        mappings = list(csv.DictReader(handle, delimiter='\t'))
    paths = [audit_path, parent_path, ledger_path, mapping_path, structure_path]
    checks = []
    for row in ledger['rows']:
        if row['status'] == 'retrieval-failed':
            checks.append({'chebi_id': row['requested_chebi_id'], 'status': 'retrieval-failed'})
            continue
        path = Path(row['path']); raw = path.read_bytes()
        if len(raw) != row['bytes'] or hashlib.sha256(raw).hexdigest() != row['sha256']:
            raise ValueError('Retrieved record checksum mismatch')
        data = json.loads(raw)
        if data['chebi_accession'] != row['returned_chebi_id']:
            raise ValueError('Returned identifier mismatch')
        paths.append(path)
        smiles = (data.get('default_structure') or {}).get('smiles')
        status = ('identifier-redirect-review' if data['chebi_accession'] != row['requested_chebi_id'] else
                  'missing-default-structure' if not smiles else
                  'invalid-or-generic-structure-excluded' if exact_key(smiles) is None else
                  'exact-id-structure-available')
        checks.append({'chebi_id': row['requested_chebi_id'], 'returned_chebi_id': data['chebi_accession'],
                       'status': status, 'source_path': str(path)})
        if status in ('exact-id-structure-available', 'invalid-or-generic-structure-excluded'):
            structures.append((data['chebi_accession'], smiles))
    report = build(json.loads(audit_path.read_text()), structures, mappings)
    report['schema'] = 'cannabis-carbon.phase1-protonation-verified-join.v1'
    report['retrieved_identity_checks'] = checks
    report['summary']['retrieved_identity_status_counts'] = dict(Counter(r['status'] for r in checks))
    report['summary']['mapping_origin_counts'] = dict(Counter(e['origin'] for r in report['rows'] for e in r['mapping_evidence']))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-protonation-verified-join.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
