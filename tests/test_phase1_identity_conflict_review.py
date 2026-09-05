import copy
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase
from cannabis_carbon.phase1_identity_conflict_review import assertions, build

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_every_conflicted_assertion_keeps_structure_and_registry_provenance(monkeypatch):
    monkeypatch.chdir(ROOT)
    audit, report = read('phase1-no-producer-audit'), read('phase1-identity-conflict-review')
    before = copy.deepcopy(audit)
    actual = build(audit, report['lookups'])
    assert actual['rdkit_version'] == rdBase.rdkitVersion
    assert {k: v for k, v in actual.items() if k != 'rdkit_version'} == {
        k: v for k, v in report.items() if k not in ('source_sha256', 'rdkit_version')}
    assert audit == before
    assert len(report['assertions']) == 142 and len(report['comparisons']) == 71
    assert report['summary']['reported_key_disagreements'] == 9
    assert report['summary']['assertions_with_terpedia_matches'] == 14
    assert report['summary']['assertions_with_pubchem_matches'] == 132
    originals = {r['id']: r for r in assertions(audit)}
    for row in report['assertions']:
        assert row['source_assertion'] == originals[row['id']]['source_assertion']
        assert row['computed_inchikey'] == Chem.MolToInchiKey(Chem.MolFromSmiles(row['source_assertion']['smiles']))
        assert all(p['InChIKey'] == row['computed_inchikey'] for p in row['pubchem_matches'])
        assert all(p['inchikey'] == row['computed_inchikey'] for p in row['terpedia_matches'])
        assert row['registry_lookup_status'] == {'terpedia': ['retrieved'], 'pubchem-keys': ['retrieved']}
    for lookup in report['lookups']:
        assert lookup['status'] == 'retrieved'
        assert hashlib.sha256((ROOT / lookup['snapshot']).read_bytes()).hexdigest() == lookup['sha256']
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_name_corroboration_does_not_prefer_one_export_globally_or_merge_tautomers():
    report = read('phase1-identity-conflict-review')
    reviews = {r['cannabisdb_id']: r for r in report['priority_reviews']}
    assert len(reviews) == 4
    assert reviews['CDB006156']['exact_structure_supported_assertion_ids'] == ['CDB006156:xml_assertion']
    assert reviews['CDB000142']['exact_structure_supported_assertion_ids'] == ['CDB000142:sdf_derived_assertion']
    assert reviews['CDB000546']['exact_structure_supported_assertion_ids'] == ['CDB000546:sdf_derived_assertion']
    assert reviews['CDB006169']['status'] == 'identity-remains-unresolved'
    rows = {r['id']: r for r in report['assertions']}
    sdf, xml = rows['CDB000546:sdf_derived_assertion'], rows['CDB000546:xml_assertion']
    assert sdf['computed_inchikey'] == xml['computed_inchikey']
    assert sdf['canonical_smiles'] != xml['canonical_smiles']
    assert sdf['pubchem_matches'][0]['exact_encoded_structure_match'] is True
    assert xml['pubchem_matches'][0]['exact_encoded_structure_match'] is False
    assert all(r['resolution_status'].startswith('unresolved;') for r in report['comparisons'])
