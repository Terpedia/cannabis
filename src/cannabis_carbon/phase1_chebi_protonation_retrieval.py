"""Bounded, resumable retrieval of exact ChEBI identifier leads (not assignments)."""
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen


def run():
    source = Path('data/reports/phase1-protonation-source-join.json')
    report = json.loads(source.read_text())
    ids = sorted({e['source_chebi_id'] for row in report['rows']
                  for e in row['unresolved_source_structure_leads']})
    if len(ids) > 500 or any(not i.startswith('CHEBI:') or not i[6:].isdigit() for i in ids):
        raise ValueError('Invalid or excessive retrieval inventory')
    folder = Path('data/raw/chebi-protonation-leads-20260906')
    folder.mkdir(exist_ok=True)

    def fetch(chebi):
        path = folder / (chebi.replace(':', '-') + '.json')
        url = 'https://www.ebi.ac.uk/chebi/backend/api/public/compound/' + chebi + '/'
        try:
            if path.exists():
                raw = path.read_bytes()
            else:
                with urlopen(Request(url, headers={'Accept': 'application/json'}), timeout=30) as response:
                    raw = response.read(5000001)
                if len(raw) > 5000000:
                    raise ValueError('Response exceeds 5 MB bound')
                data = json.loads(raw)
                if not isinstance(data, dict) or 'chebi_accession' not in data:
                    raise ValueError('Unexpected response')
                with path.open('xb') as handle:
                    handle.write(raw)
            data = json.loads(raw)
            return {'requested_chebi_id': chebi, 'returned_chebi_id': data['chebi_accession'],
                    'url': url, 'path': str(path), 'sha256': hashlib.sha256(raw).hexdigest(),
                    'bytes': len(raw), 'has_default_smiles': bool((data.get('default_structure') or {}).get('smiles')),
                    'status': 'retrieved-exact-id' if data['chebi_accession'] == chebi else 'identifier-redirect-review'}
        except Exception as exc:
            return {'requested_chebi_id': chebi, 'url': url, 'status': 'retrieval-failed', 'error': str(exc)}

    rows = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        for row in pool.map(fetch, ids):
            rows.append(row)
            if len(rows) % 50 == 0:
                print(f'ChEBI retrieval: {len(rows)}/{len(ids)}', flush=True)
    output = {'retrieved_at': datetime.now(timezone.utc).isoformat(),
              'source_sha256': {str(source): hashlib.sha256(source.read_bytes()).hexdigest()},
              'rows': rows, 'claim_boundary': 'Identifier leads only; returned structures require exact '
              'comparison before target assignment. Redirects, missing structures and failures remain explicit.'}
    (folder / 'retrievals.json').write_text(json.dumps(output, separators=(',', ':')) + '\n')
    print(json.dumps({'requested': len(rows), 'failed': sum(r['status'] == 'retrieval-failed' for r in rows)}), flush=True)


if __name__ == '__main__':
    run()
