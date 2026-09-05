"""Read-only verification of the immutable Phase 1 thiolase publication batch."""
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path

NAMES = ('synthase-precursor-audit', 'cbga-input-audit', 'evidence-weighted-routes',
         'weighted-gap-search', 'thiolase-candidate-net')
PROJECT = 'terpedia-489015'
DATASET = 'terpedia_core'
PRINCIPAL = 'cannabis-metabolome@terpedia-489015.iam.gserviceaccount.com'


def bq(*args):
    result = subprocess.run(['/Users/danielmcshan/google-cloud-sdk/bin/bq',
        '--project_id=' + PROJECT, '--location=us-central1', '--format=json', *args],
        capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def canonical(rows):
    return Counter(json.dumps(r, sort_keys=True, separators=(',', ':')) for r in rows)


def verify():
    tables = []
    for name in NAMES:
        suffix = name.replace('-', '_') + '_20260905_v1'
        table = 'cannabis_phase1_' + suffix
        job_id = 'cannabis_sa_' + suffix
        path = Path('data/reports/phase1-' + name + '.json')
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        rows = [json.loads(line) for line in Path('data/derived/phase1-' + name + '.ndjson').read_text().splitlines()]
        if len({(r['record_kind'], r['record_id']) for r in rows}) != len(rows):
            raise ValueError('Duplicate local record identity')
        if {r['report_sha256'] for r in rows} != {sha}:
            raise ValueError('Local export does not match report hash')
        job = bq('show', '--job=true', job_id)
        status = job['status']
        destination = job['configuration']['load']['destinationTable']
        if status['state'] != 'DONE' or status.get('errorResult') or status.get('errors'):
            raise ValueError('Load not successfully completed')
        if job['user_email'] != PRINCIPAL or destination != {
                'projectId': PROJECT, 'datasetId': DATASET, 'tableId': table}:
            raise ValueError('Unexpected principal or destination')
        remote = bq('query', '--use_legacy_sql=false', '--max_rows=20000',
            '--maximum_bytes_billed=33554432', 'SELECT * FROM `' + PROJECT + '.' + DATASET + '.' + table + '`')
        if canonical(rows) != canonical(remote):
            raise ValueError('Full stored record multiset differs from local export')
        tables.append({'table': table, 'job_id': job_id, 'rows': len(rows), 'report_sha256': sha})
        print('Verified ' + table + ': ' + str(len(rows)) + ' rows', flush=True)
    return {'schema': 'cannabis-carbon.gcp-batch-verification.v1', 'project': PROJECT,
        'dataset': DATASET, 'principal': PRINCIPAL, 'verified': True,
        'verification': 'All load jobs DONE without errors; exact service-account principal and destination checked. Full bounded SELECT readbacks match complete local NDJSON record multisets, report checksums and unique kind/id pairs.',
        'tables': tables}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
