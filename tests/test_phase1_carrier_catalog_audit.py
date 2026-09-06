import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_carrier_catalog_audit import index_sources, context_key
from cannabis_carbon.phase1_glycerophospholipid_input_audit import assemble_current


def test_full_carrier_catalog_joins_and_inventory_impact():
    root = Path('data/reports'); raw = Path('data/raw/rhea-carrier-snapshot-20260906')
    read = lambda n: json.loads((root / ('phase1-' + n + '.json')).read_bytes())
    audit = read('carrier-catalog-audit'); current = read('geranial-reduction-net')
    for p, sha in audit['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    properties = json.loads((raw / 'carrier-properties.json').read_bytes())['results']['bindings']
    occurrences = json.loads((raw / 'carrier-occurrences.json').read_bytes())['results']['bindings']
    carriers, source_index = index_sources(properties, occurrences)
    assert audit['carriers'] == carriers and len(carriers) == 1980
    assert len(source_index) == 6730 and sum(map(len, source_index.values())) == 13974
    layers = [json.loads((root / n).read_bytes()) for n in current['baseline_certificate_reports']]
    reactions, _, _ = assemble_current(read('full-balanced-network'), read('marts-completions'), current, layers)
    expected = {rid for rid, r in reactions.items() if any(s.get('source_reaction_id') in source_index for s in r.get('sources', []))}
    assert {r['id'] for r in audit['equations']} == expected
    assert len(expected) == len(audit['equations']) == 1316
    for row in audit['equations']:
        r = reactions[row['id']]; assert row['model_reaction'] == r
        joins = [{'source_record': s, 'carrier_occurrences': source_index[s['source_reaction_id']]}
                 for s in r.get('sources', []) if s.get('source_reaction_id') in source_index]
        assert row['carrier_source_joins'] == joins
        assert row['distinct_unoriented_carrier_contexts'] == len({context_key(j['carrier_occurrences']) for j in joins})
        assert row['requires_carrier_specific_reconstruction']
    certs = {c['compound_id']: c for layer in [*layers, current] for c in layer.get('certificates', []) + layer.get('new_certificates', [])}
    affected = set()
    assert len(audit['targets']) == len(current['targets']) == 6220
    for row, target in zip(audit['targets'], current['targets']):
        for k in ('cannabisdb_id', 'label', 'compound_id', 'net_status'):
            assert row[k] == target[k]
        cid = row['compound_id']; hits = sorted({s['reaction_id'] for s in certs.get(cid, {}).get('steps', [])} & expected)
        assert row['carrier_source_equation_ids'] == hits
        if hits:
            affected.add(cid)
            assert row['carrier_audit_status'] == 'carrier-source-reconstruction-required'
        else:
            assert 'not-full-identity-validation' in row['carrier_audit_status'] or cid not in certs
    assert len(affected) == audit['summary']['affected_certificates'] == 2669
    assert sum(bool(t['carrier_source_equation_ids']) for t in audit['targets']) == audit['summary']['affected_inventory_records'] == 2672
    assert sum(r['distinct_unoriented_carrier_contexts'] > 1 for r in audit['equations']) == audit['summary']['model_equations_merging_multiple_carrier_contexts'] == 10
    assert audit['summary']['model_equations'] == len(reactions) == 18222
    assert audit['summary']['screened_certificates'] == len(certs) == 2737
    assert not audit['summary']['corrected_pathway_coverage_established']
