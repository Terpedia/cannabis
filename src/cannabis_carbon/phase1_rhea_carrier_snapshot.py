"""Retrieve full Rhea carrier identities and their directed-reaction occurrences."""
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ENDPOINT = 'https://sparql.rhea-db.org/sparql'
PREFIX = 'PREFIX rh: <http://rdf.rhea-db.org/> '
QUERIES = {
    'carrier-properties': PREFIX + '''SELECT DISTINCT ?subject ?predicate ?object WHERE {
      { ?subject rh:reactivePart ?part . ?subject ?predicate ?object }
      UNION
      { ?carrier rh:reactivePart ?subject . ?subject ?predicate ?object }
    } ORDER BY ?subject ?predicate ?object''',
    'carrier-occurrences': PREFIX + '''SELECT DISTINCT ?reaction ?sideRole ?side ?participant ?carrier ?coefficientPredicate WHERE {
      VALUES ?sideRole { rh:substrates rh:products }
      ?reaction ?sideRole ?side .
      ?side rh:contains ?participant ; ?coefficientPredicate ?participant .
      ?participant rh:compound ?carrier .
      ?carrier rh:reactivePart ?part .
      FILTER(?coefficientPredicate != rh:contains)
    } ORDER BY ?reaction ?sideRole ?side ?participant ?carrier ?coefficientPredicate'''
}


def run():
    root = Path('data/raw/rhea-carrier-snapshot-20260906')
    root.mkdir(exist_ok=True)
    receipts = []
    for name, query in QUERIES.items():
        path = root / (name + '.json')
        if path.exists():
            raise ValueError('Snapshot exists; inspect it rather than overwrite or refetch')
        url = ENDPOINT + '?' + urllib.parse.urlencode({'query': query, 'format': 'json'})
        request = urllib.request.Request(url, headers={'Accept': 'application/sparql-results+json',
            'User-Agent': 'Terpedia-Cannabis-carrier-identity-audit/1.0'})
        with urllib.request.urlopen(request, timeout=55) as response:
            payload = response.read()
            headers = dict(response.headers)
        result = json.loads(payload)
        rows = result['results']['bindings']
        if not rows:
            raise ValueError('Empty carrier snapshot')
        if name == 'carrier-properties' and not any(
                r['subject']['value'] == 'http://rdf.rhea-db.org/Compound_12863'
                and r['predicate']['value'] == 'http://rdf.rhea-db.org/reactivePart' for r in rows):
            raise ValueError('Known cytochrome carrier missing')
        if name == 'carrier-occurrences' and not any(
                r['reaction']['value'] == 'http://rdf.rhea-db.org/56273'
                and r['carrier']['value'] == 'http://rdf.rhea-db.org/Compound_12863' for r in rows):
            raise ValueError('Known cytochrome reaction occurrence missing')
        path.write_bytes(payload)
        receipt = {'name': name, 'endpoint': ENDPOINT, 'query': query,
            'retrieved_at': datetime.now(timezone.utc).isoformat(), 'headers': headers,
            'path': str(path), 'bytes': len(payload), 'rows': len(rows),
            'sha256': hashlib.sha256(payload).hexdigest()}
        (root / (name + '-retrieval.json')).write_text(json.dumps(receipt, indent=2) + '\n')
        receipts.append(receipt)
        print(json.dumps({k: receipt[k] for k in ('name', 'rows', 'bytes', 'sha256')}), flush=True)
    (root / 'manifest.json').write_text(json.dumps({'retrievals': receipts,
        'claim_boundary': 'Current endpoint snapshot of reactive-part carriers, not a complete polymer audit or a corrected metabolic model. Retain full source carrier URIs and reactive-part identities separately. Completeness must be checked against independent endpoint counts and the pinned historical reaction inventory before use.'}, indent=2) + '\n')


if __name__ == '__main__':
    run()
