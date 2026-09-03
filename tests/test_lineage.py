import gzip
import json

from cannabis_carbon.lineage import build_carbon_lineage


def test_lineage_marks_unmatched_cannabisdb_compounds(tmp_path):
    network = {"entities": [{"id": "chebi:16526", "type": "metabolite", "attributes": {"canonicalSmiles": "O=C=O"}}], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle:
        json.dump(network, handle)
    mapping_path = tmp_path / "mapping.json"
    mapping_path.write_text(json.dumps({"reactions": []}))
    crosswalk_path = tmp_path / "crosswalk.json"
    crosswalk_path.write_text(json.dumps({"matches": []}))
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": [{"id": "CDB1", "carbon_atom_count": 2}]}))
    result = build_carbon_lineage(network_path, mapping_path, crosswalk_path, compounds_path, tmp_path / "out.json")
    assert result["target_summary"] == {"supported": 0, "candidate": 0, "unresolved": 1}


def test_lineage_accepts_direction_override_file(tmp_path):
    directions = tmp_path / "directions.json"
    directions.write_text(json.dumps({"rhea:test": {"orientation": "reverse_master", "directional_rhea_id": "rhea:test-forward", "source": "test"}}))
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump({"entities": [], "statements": []}, handle)
    for name, value in (("mapping.json", {"reactions": []}), ("crosswalk.json", {"matches": []}), ("compounds.json", {"compounds": []})):
        (tmp_path / name).write_text(json.dumps(value))
    output = tmp_path / "out.json"
    build_carbon_lineage(network_path, tmp_path / "mapping.json", tmp_path / "crosswalk.json", tmp_path / "compounds.json", output, directions)
    assert json.loads(output.read_text())["direction_overrides"]["rhea:test"]["directional_rhea_id"] == "rhea:test-forward"


def test_lineage_uses_rdkit_carbon_atom_indices_for_targets(tmp_path):
    co2 = {"id": "chebi:16526", "type": "metabolite", "attributes": {"canonicalSmiles": "O=C=O"}}
    network = {
        "entities": [co2, {"id": "rhea:test", "type": "biochemical_reaction"}],
        "statements": [
            {"subjectId": "rhea:test", "predicate": "has_reactant", "objectEntityId": "chebi:16526"},
            {"subjectId": "rhea:test", "predicate": "has_product", "objectEntityId": "chebi:16526"},
        ],
    }
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    (tmp_path / "mapping.json").write_text(json.dumps({"reactions": [{"reaction_id": "rhea:test", "reaction_smiles": "O=C=O>>O=C=O", "mappings": [{"status": "inferred", "reactant_index": 0, "reactant_atom": 1, "product_index": 0, "product_atom": 1}]}]}))
    (tmp_path / "crosswalk.json").write_text(json.dumps({"matches": [{"cannabisdb": {"cannabisdb_id": "CDB1"}, "terpedia_id": "chebi:16526"}]}))
    (tmp_path / "compounds.json").write_text(json.dumps({"compounds": [{"id": "CDB1", "carbon_atom_count": 1}]}))
    result = build_carbon_lineage(network_path, tmp_path / "mapping.json", tmp_path / "crosswalk.json", tmp_path / "compounds.json", tmp_path / "out.json")
    assert result["target_summary"] == {"supported": 1, "candidate": 0, "unresolved": 0}
