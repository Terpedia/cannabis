import gzip
import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_pg_source_identity_audit import build


def test_original_xml_assertions_are_preserved_for_all_affected_accessions():
    path = Path('data/reports/phase1-pg-source-identity-audit.json')
    r = json.loads(path.read_bytes())
    for p,sha in r['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    review = json.loads(Path('data/reports/phase1-acyl-input-identity-review.json').read_bytes())
    current = json.loads(Path('data/reports/phase1-glycerophospholipid-net.json').read_bytes())
    xml = gzip.decompress(Path('data/terpedia/cannabisdb-compounds.xml.gz').read_bytes()).decode()
    assert build(review,current,xml) == {k:v for k,v in r.items() if k!='source_sha256'}
    assert r['summary']['source_records'] == 23
    assert r['summary']['new_reactions'] == r['summary']['new_CO2_route_claims'] == 0
    for t in r['targets']:
        assert t['source_fields']['accession'] == t['cannabisdb_id']
        assert t['source_name_geometry_tokens']
        assert t['unassigned_smiles_alkene_bonds'] > 0
        assert t['xml_smiles_matches_retained_structure']
        assert t['name_structure_status'] == 'name-specifies-geometry-omitted-from-structure'
        assert 'expected to be in Cannabis' in t['source_fields']['description']
        assert t['occurrence_assertion'] == 'source-says-expected-in-Cannabis-not-detection'
