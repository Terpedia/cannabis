"""Registry corroboration of disputed source assertions; never automatic repair."""
import hashlib
import json
import os
import subprocess
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors

from .phase1_target_coverage import encoded_structure
from .phase1_scope import write_rows
from .pubchem import PROPERTIES

RAW = Path('data/raw/phase1-identity-conflict-review')
TABLE = 'terpedia-489015.terpedia_core.terpene_identity_set'


def assertions(audit):
    rows = []
    for conflict in audit['source_identity_conflicts']:
        for kind in ('sdf_derived_assertion', 'xml_assertion'):
            source = conflict[kind]
            molecule = Chem.MolFromSmiles(source['smiles'])
            if molecule is None:
                raise ValueError('Invalid conflict structure')
            key = Chem.MolToInchiKey(molecule)
            if not key:
                raise ValueError('Conflict structure lacks computable InChIKey')
            rows.append({'id': conflict['cannabisdb_id'] + ':' + kind,
                'cannabisdb_id': conflict['cannabisdb_id'], 'source_kind': kind, 'source_assertion': source,
                'canonical_smiles': encoded_structure(source['smiles'])[0],
                'computed_inchikey': key, 'computed_formula': rdMolDescriptors.CalcMolFormula(molecule),
                'computed_carbon_count': sum(a.GetAtomicNum() == 6 for a in molecule.GetAtoms()),
                'reported_inchikey_matches_structure': source.get('inchikey') == key,
                'reported_formula_matches_computed_text': source.get('formula') == rdMolDescriptors.CalcMolFormula(molecule),
                'xref_origin': 'CannabisDB XML enrichment; the two copies are not independent corroboration'})
    return rows


def fetch(kind, request, keys=None):
    RAW.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256((kind + request).encode()).hexdigest()
    cache = RAW / (digest + '.request.json')
    if cache.exists():
        record = json.loads(cache.read_text())
        if record['request'] != request or record['kind'] != kind:
            raise ValueError('Request cache mismatch')
        if record['status'] == 'retrieved' and hashlib.sha256(Path(record['snapshot']).read_bytes()).hexdigest() != record['sha256']:
            raise ValueError('Response checksum mismatch')
        return record
    record = {'kind': kind, 'request': request, 'requested_keys': keys or [],
              'retrieved_at': datetime.now(timezone.utc).isoformat()}
    try:
        if kind == 'terpedia':
            command = [os.environ.get('BQ_BINARY', 'bq'), '--project_id=terpedia-489015', '--location=us-central1',
                       '--format=prettyjson', 'query', '--use_legacy_sql=false', '--maximum_bytes_billed=536870912',
                       '--max_rows=10000', request]
            payload = subprocess.check_output(command)
        else:
            with urllib.request.urlopen(request, timeout=45) as response:
                payload = response.read()
        json.loads(payload)
        snapshot = RAW / (digest + '.response.json')
        snapshot.write_bytes(payload)
        record.update(status='retrieved', snapshot=str(snapshot), sha256=hashlib.sha256(payload).hexdigest())
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        record.update(status='not-found' if isinstance(error, urllib.error.HTTPError) and error.code == 404 else 'retrieval-failed', reason=str(error))
    cache.write_text(json.dumps(record, separators=(',', ':')) + '\n')
    return record


def build(audit, lookups):
    rows = assertions(audit)
    terpedia, pubchem, name_records = [], [], {}
    for lookup in lookups:
        if lookup['status'] != 'retrieved':
            continue
        payload = Path(lookup['snapshot']).read_bytes()
        if hashlib.sha256(payload).hexdigest() != lookup['sha256']:
            raise ValueError('Registry response checksum mismatch')
        data = json.loads(payload)
        if lookup['kind'] == 'terpedia':
            if any(r['inchikey'] not in lookup['requested_keys'] for r in data):
                raise ValueError('Unrequested Terpedia identity')
            terpedia.extend(data)
        elif lookup['kind'] == 'pubchem-keys':
            returned = data['PropertyTable']['Properties']
            if any(r['InChIKey'] not in lookup['requested_keys'] for r in returned):
                raise ValueError('Unrequested PubChem identity')
            pubchem.extend(returned)
        elif lookup['kind'] == 'pubchem-name':
            name_records[lookup['requested_keys'][0]] = data['PropertyTable']['Properties']
    for row in rows:
        key = row['computed_inchikey']
        row['registry_lookup_status'] = {kind: [l['status'] for l in lookups
            if l['kind'] == kind and key in l['requested_keys']] for kind in ('terpedia', 'pubchem-keys')}
        row['terpedia_matches'] = [r for r in terpedia if r['inchikey'] == key]
        row['pubchem_matches'] = []
        for p in pubchem:
            if p['InChIKey'] != key:
                continue
            smiles = p.get('SMILES') or p.get('IsomericSMILES')
            row['pubchem_matches'].append({**p, 'exact_encoded_structure_match':
                encoded_structure(smiles)[0] == row['canonical_smiles'] if smiles else None,
                'source_url': 'https://pubchem.ncbi.nlm.nih.gov/compound/' + str(p['CID'])})
        name = row['source_assertion'].get('name') or row['source_assertion'].get('label')
        row['priority_name_lookup_matches'] = name_records.get(name, [])
        row['name_query_exact_structure_cids'] = sorted({p['CID'] for p in name_records.get(name, [])
            if encoded_structure(p.get('SMILES') or p.get('IsomericSMILES'))[0] == row['canonical_smiles']})
        row['status'] = 'source-assertion-under-review; registry-match-not-correction'
    by_id = {}
    for row in rows:
        by_id.setdefault(row['cannabisdb_id'], {})[row['source_kind']] = row
    comparisons = []
    priority_reviews = []
    for cid, pair in by_id.items():
        a, b = pair['sdf_derived_assertion'], pair['xml_assertion']
        ka, kb = a['computed_inchikey'], b['computed_inchikey']
        relationship = 'same-standard-InChIKey; encoded-structures-still-distinct' if ka == kb else 'same-connectivity-key-only' if ka.split('-')[0] == kb.split('-')[0] else 'different-connectivity-key'
        comparisons.append({'cannabisdb_id': cid, 'assertion_ids': [a['id'], b['id']], 'relationship': relationship,
                            'resolution_status': 'unresolved; no source selected or merged'})
        name = a['source_assertion'].get('label')
        if name in name_records:
            choices = [r for r in (a, b) if r['name_query_exact_structure_cids']]
            priority_reviews.append({'cannabisdb_id': cid, 'queried_name': name,
                'exact_structure_supported_assertion_ids': [r['id'] for r in choices],
                'pubchem_name_query_records': name_records[name],
                'status': 'provisional-registry-supported-source-choice-for-name' if len(choices) == 1 else 'identity-remains-unresolved',
                'action': 'Prepare explicit versioned source reconciliation; do not overwrite historical structures or xrefs.' if len(choices) == 1 else 'Review stereochemistry and source primary records before selecting a structure.',
                'claim_boundary': 'Registry exact-structure/name corroboration, not Cannabis occurrence, primary source correction or automatic migration.'})
    return {'schema': 'cannabis-carbon.phase1-identity-conflict-review.v1', 'assertions': rows,
            'comparisons': comparisons, 'priority_reviews': priority_reviews, 'lookups': lookups,
            'summary': {'conflicted_accessions': len(comparisons), 'source_assertions': len(rows),
                        'relationship_counts': dict(Counter(r['relationship'] for r in comparisons)),
                        'assertions_with_terpedia_matches': sum(bool(r['terpedia_matches']) for r in rows),
                        'assertions_with_pubchem_matches': sum(bool(r['pubchem_matches']) for r in rows),
                        'reported_key_disagreements': sum(not r['reported_inchikey_matches_structure'] for r in rows),
                        'failed_lookups': sum(l['status'] == 'retrieval-failed' for l in lookups),
                        'priority_name_reviews': len(priority_reviews),
                        'provisional_named_structure_choices': sum(r['status'] == 'provisional-registry-supported-source-choice-for-name' for r in priority_reviews)},
            'rdkit_version': rdBase.rdkitVersion,
            'claim_boundary': 'Registry corroboration of both source structures, not source selection, confirmed names, Cannabis occurrence or pathways. Standard InChIKey equality can normalize distinctions; exact encoded SMILES comparisons remain separate. Name lookup is an independent query, not authoritative structure assignment. Source xrefs are shared XML assertions, not independent SDF/XML agreement. No graph identity or carbon accounting changed. Atom tracing remains deferred.'}


def run():
    RDLogger.DisableLog('rdApp.warning')
    path = Path('data/reports/phase1-no-producer-audit.json')
    audit = json.loads(path.read_text())
    for source, sha in audit['source_sha256'].items():
        if hashlib.sha256(Path(source).read_bytes()).hexdigest() != sha:
            raise ValueError('Source lineage changed')
    rows = assertions(audit)
    keys = sorted({r['computed_inchikey'] for r in rows})
    sql = ('SELECT terpene_id,inchikey,smiles,molecular_formula,source_memberships,source_record_ids,source_crossrefs '
           'FROM `' + TABLE + '` WHERE inchikey IN (' + ','.join("'" + k + "'" for k in keys) + ') ORDER BY inchikey,terpene_id')
    lookups = [fetch('terpedia', sql, keys)]
    for start in range(0, len(keys), 25):
        batch = keys[start:start + 25]
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/' + ','.join(batch) + '/property/' + PROPERTIES + '/JSON'
        lookups.append(fetch('pubchem-keys', url, batch))
    for name in ('Ribitol', 'D-arabitol', 'Acetamide', 'Glycerol'):
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' + urllib.parse.quote(name, safe='') + '/property/' + PROPERTIES + '/JSON'
        lookups.append(fetch('pubchem-name', url, [name]))
    report = build(audit, lookups)
    report['source_sha256'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-identity-conflict-review.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('assertion', 'assertions', 'id'), ('comparison', 'comparisons', 'cannabisdb_id'), ('priority_review', 'priority_reviews', 'cannabisdb_id'), ('lookup', 'lookups', 'request')]
    records = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, collection, identifier in groups:
        records.extend((kind, r[identifier], r) for r in report[collection])
    count = write_rows(records, sha, Path('data/derived/phase1-identity-conflict-review.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
