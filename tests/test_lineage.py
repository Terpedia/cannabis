import gzip
import json

from rdkit import Chem

from cannabis_carbon.lineage import _entity_atom_index_map, build_carbon_atom_audit, build_carbon_lineage


def test_lineage_remaps_reaction_atom_order_to_entity_atom_order():
    reaction_molecule = Chem.MolFromSmiles("O=C(O)C")
    atom_map = _entity_atom_index_map(reaction_molecule, "CC(=O)O")
    assert atom_map is not None
    assert atom_map[1] == 1
    assert atom_map[3] == 0


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


def test_direction_evidence_retains_uniprot_directional_rhea_id():
    from cannabis_carbon.lineage import _direction_evidence

    network = {"statements": [{
        "subjectId": "rhea:test", "predicate": "physiological_direction_right_to_left",
        "sources": [{"url": "https://example.test/uniprot"}],
        "qualifiers": {"support": [{"directionalRheaId": "12346"}]},
    }]}
    directions, conflicts = _direction_evidence(network)
    assert not conflicts
    assert directions["rhea:test"] == {"directional_rhea_id": "12346", "orientation": "reverse_master", "source": "https://example.test/uniprot", "reason": "Terpedia physiological direction statement"}


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


def test_lineage_marks_reachable_connectivity_identity_as_candidate(tmp_path):
    network = {"entities": [{"id": "chebi:16526", "type": "metabolite", "attributes": {"canonicalSmiles": "O=C=O"}}], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    for name, value in (("mapping.json", {"reactions": []}), ("crosswalk.json", {"matches": [], "candidate_matches": [{"terpedia_id": "chebi:16526", "terpedia_label": "CO2", "method": "connectivity-inchikey-candidate", "cannabisdb": {"cannabisdb_id": "CDB1"}}]}), ("compounds.json", {"compounds": [{"id": "CDB1", "carbon_atom_count": 1}]})):
        (tmp_path / name).write_text(json.dumps(value))
    result = build_carbon_lineage(network_path, tmp_path / "mapping.json", tmp_path / "crosswalk.json", tmp_path / "compounds.json", tmp_path / "out.json")
    assert result["target_summary"] == {"supported": 0, "candidate": 1, "unresolved": 0}


def test_carbon_atom_audit_accounts_for_every_target_carbon(tmp_path):
    network = {"entities": [{"id": "chebi:16526", "type": "metabolite", "attributes": {"canonicalSmiles": "O=C=O"}}], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle:
        json.dump(network, handle)
    lineage_path = tmp_path / "lineage.json"
    lineage_path.write_text(json.dumps({"co2_entity_id": "chebi:16526", "carbon_source_policy": "CO2 only", "carbon_edges": [], "targets": [{"cannabisdb_id": "CDB1", "terpedia_id": "chebi:16526", "identity_status": "exact", "status": "supported", "reason": "co2-seed"}]}))
    crosswalk_path = tmp_path / "crosswalk.json"
    crosswalk_path.write_text(json.dumps({"matches": [{"cannabisdb": {"cannabisdb_id": "CDB1"}, "terpedia_id": "chebi:16526"}]}))
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": [{"id": "CDB1", "smiles": "O=C=O", "carbon_atom_count": 1, "source_url": "test"}]}))
    output = tmp_path / "audit.json"
    result = build_carbon_atom_audit(network_path, lineage_path, crosswalk_path, compounds_path, output)
    assert result["carbon_atoms_total"] == 1
    assert result["status_counts"] == {"supported": 1, "candidate": 0, "inferred": 0, "unresolved": 0}
    groups = json.loads(output.read_text())["compounds"][0]["groups"]
    assert sum(len(group["atom_indices"]) for group in groups) == 1
    assert json.loads(output.read_text())["compounds"][0]["co2_paths"] == {"1": []}
