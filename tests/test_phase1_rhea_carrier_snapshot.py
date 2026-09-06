import hashlib
import json
from pathlib import Path


def test_source_carriers_occurrences_counts_and_identity_links():
    root = Path('data/raw/rhea-carrier-snapshot-20260906')
    read = lambda name: json.loads((root / name).read_bytes())
    manifest = read('manifest.json')
    for receipt in manifest['retrievals']:
        payload = Path(receipt['path']).read_bytes()
        assert len(payload) == receipt['bytes']
        assert hashlib.sha256(payload).hexdigest() == receipt['sha256']
        assert len(json.loads(payload)['results']['bindings']) == receipt['rows']
        assert receipt == read(receipt['name'] + '-retrieval.json')
    rh = 'http://rdf.rhea-db.org/'
    properties = read('carrier-properties.json')['results']['bindings']
    occurrences = read('carrier-occurrences.json')['results']['bindings']
    carriers = {r['subject']['value'] for r in properties if r['predicate']['value'] == rh + 'reactivePart'}
    counts = {r['name']: r['result'] for r in read('independent-counts.json')['counts']}
    assert len(carriers) == counts['distinct-carriers-with-reactive-parts'] == 1980
    assert len(occurrences) == counts['directed-carrier-occurrence-rows'] == 13974
    assert {r['carrier']['value'] for r in occurrences} == carriers
    assert len({json.dumps(r, sort_keys=True) for r in occurrences}) == len(occurrences)
    assert all(r['sideRole']['value'] in (rh + 'substrates', rh + 'products') for r in occurrences)
    # Variable stoichiometry is source data, not an integer coefficient.
    assert any(r['coefficientPredicate']['value'] == rh + 'contains2n' for r in occurrences)
    for carrier in ('12863', '12864', '14399', '10350'):
        assert rh + 'Compound_' + carrier in carriers
    assert any(r['reaction']['value'] == rh + '56273' and r['carrier']['value'] == rh + 'Compound_12863' for r in occurrences)
    assert any(r['reaction']['value'] == rh + '79540' and r['carrier']['value'] == rh + 'Compound_14399' for r in occurrences)
