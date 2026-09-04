import json

from cannabis_carbon.validate import validate_artifacts


def test_validate_artifacts_enforces_carbon_partition_and_mapping_rows(tmp_path):
    compounds = tmp_path / "compounds.json"
    compounds.write_text(json.dumps({"compounds": [{"id": "CDB1", "smiles": "CC"}]}))
    audit = tmp_path / "audit.json"
    audit.write_text(json.dumps({"carbon_atoms_total": 2, "compounds": [{"cannabisdb_id": "CDB1", "groups": [{"status": "unresolved", "reason": "test", "provenance": ["test"] , "atom_indices": [0]}, {"status": "unresolved", "reason": "test", "provenance": ["test"], "atom_indices": [1]}]}]}))
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"reactions": [{"reaction_id": "r1", "product_carbon_atom_count": 1, "mappings": [{"product_atom": 0, "status": "inferred"}]}]}))
    balance = tmp_path / "balance.json"
    balance.write_text(json.dumps({"reactions": [{"reaction_id": "r1", "status": "balanced"}]}))
    result = validate_artifacts(audit, mapping, balance, compounds, tmp_path / "out.json")
    assert result["valid"]
    assert result["checks"] == {"carbon_atom_partition": True, "atom_evidence_fields": True, "global_carbon_accounting": True, "carbon_path_integrity": True, "carbon_mapping_queue_consistency": True, "reaction_product_carbon_rows": True, "reaction_mapping_classification": True, "no_imbalanced_reactions": True}


def test_validate_artifacts_rejects_omitted_carbon(tmp_path):
    compounds = tmp_path / "compounds.json"
    compounds.write_text(json.dumps({"compounds": [{"id": "CDB1", "smiles": "CC"}]}))
    audit = tmp_path / "audit.json"
    audit.write_text(json.dumps({"carbon_atoms_total": 1, "compounds": [{"cannabisdb_id": "CDB1", "groups": [{"atom_indices": [0]}]}]}))
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"reactions": []}))
    balance = tmp_path / "balance.json"
    balance.write_text(json.dumps({"reactions": []}))
    result = validate_artifacts(audit, mapping, balance, compounds, tmp_path / "out.json")
    assert not result["valid"]
    assert result["failures"][0]["kind"] == "carbon-atom-partition-mismatch"


def test_validate_artifacts_rejects_atom_without_evidence_fields(tmp_path):
    compounds = tmp_path / "compounds.json"
    compounds.write_text(json.dumps({"compounds": [{"id": "CDB1", "smiles": "CC"}]}))
    audit = tmp_path / "audit.json"
    audit.write_text(json.dumps({"carbon_atoms_total": 2, "compounds": [{"cannabisdb_id": "CDB1", "groups": [{"status": "unresolved", "provenance": ["test"], "atom_indices": [0, 1]}]}]}))
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"reactions": []}))
    balance = tmp_path / "balance.json"
    balance.write_text(json.dumps({"reactions": []}))
    result = validate_artifacts(audit, mapping, balance, compounds, tmp_path / "out.json")
    assert not result["valid"]
    assert any(f["kind"] == "atom-evidence-fields-missing" for f in result["failures"])
