import gzip
import json

from cannabis_carbon.networkdb import build_map_snapshot, build_networkdb


def test_networkdb_contains_all_source_records(tmp_path):
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle:
        json.dump({"entities": [{"id": "m:1", "type": "metabolite", "label": "M"}, {"id": "r:1", "type": "biochemical_reaction", "attributes": {}, "identifiers": {}}], "statements": []}, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": [{"id": "CDB1", "label": "C"}]}))
    crosswalk_path = tmp_path / "crosswalk.json"
    crosswalk_path.write_text(json.dumps({"matches": [], "ambiguous": 0, "unmatched": 1}))
    output = tmp_path / "networkdb.json"
    coverage = build_networkdb(network_path, compounds_path, crosswalk_path, output)
    result = json.loads(output.read_text())
    assert coverage["cannabisdb_compounds"] == 1
    assert coverage["terpedia_metabolites"] == 1
    assert coverage["terpedia_reactions"] == 1
    assert len(result["compounds"]) == 2
    assert len(result["reactions"]) == 1
    assert result["reactions"][0]["carbon_mapping"]["status"] == "unavailable"
    assert result["reactions"][0]["carbon_mapping"]["lineage_edge_counts"] == {"inferred": 0, "candidate": 0}


def test_networkdb_reads_object_reaction_enzyme_associations(tmp_path):
    network_path = tmp_path / "network.json.gz"
    network = {"entities": [
        {"id": "r:1", "type": "biochemical_reaction", "label": "A = B", "attributes": {}, "identifiers": {}},
        {"id": "p:1", "type": "protein", "label": "enzyme"},
    ], "statements": [{
        "subjectId": "p:1", "predicate": "has_catalytic_activity", "objectEntityId": "r:1",
        "qualifiers": {"directExperimentalEvidence": False}, "sources": [{"url": "https://example.test"}],
    }]}
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    crosswalk_path = tmp_path / "crosswalk.json"
    crosswalk_path.write_text(json.dumps({"matches": [], "ambiguous": 0, "unmatched": 0}))
    output = tmp_path / "networkdb.json"
    build_networkdb(network_path, compounds_path, crosswalk_path, output)
    reaction = json.loads(output.read_text())["reactions"][0]
    assert reaction["enzyme_ids"] == ["p:1"]
    assert reaction["enzyme_associations"][0]["predicate"] == "has_catalytic_activity"
    assert reaction["status"] == "candidate"


def test_networkdb_labels_non_enzymatic_reactions(tmp_path):
    network_path = tmp_path / "network.json.gz"
    network = {"entities": [{
        "id": "r:1", "type": "biochemical_reaction", "label": "decarboxylation",
        "attributes": {"reactionClass": "non-enzymatic-decarboxylation"}, "identifiers": {},
    }], "statements": []}
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    crosswalk_path = tmp_path / "crosswalk.json"
    crosswalk_path.write_text(json.dumps({"matches": [], "ambiguous": 0, "unmatched": 0}))
    output = tmp_path / "networkdb.json"
    build_networkdb(network_path, compounds_path, crosswalk_path, output)
    assert json.loads(output.read_text())["reactions"][0]["status"] == "non_enzymatic"


def test_networkdb_orients_participants_from_directional_override(tmp_path):
    network_path = tmp_path / "network.json.gz"
    network = {"entities": [
        {"id": "m:a", "type": "metabolite", "label": "A"},
        {"id": "m:b", "type": "metabolite", "label": "B"},
        {"id": "r:1", "type": "biochemical_reaction", "label": "B = A", "attributes": {}, "identifiers": {}},
    ], "statements": [
        {"subjectId": "r:1", "predicate": "has_reactant", "objectEntityId": "m:b"},
        {"subjectId": "r:1", "predicate": "has_product", "objectEntityId": "m:a"},
    ]}
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    (tmp_path / "directional-reaction-overrides.json").write_text(json.dumps({"r:1": {"orientation": "reverse_master", "directional_rhea_id": "2"}}))
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    crosswalk_path = tmp_path / "crosswalk.json"
    crosswalk_path.write_text(json.dumps({"matches": [], "ambiguous": 0, "unmatched": 0}))
    output = tmp_path / "networkdb.json"
    build_networkdb(network_path, compounds_path, crosswalk_path, output)
    reaction = json.loads(output.read_text())["reactions"][0]
    assert [p["compound_id"] for p in reaction["reactants"]] == ["m:a"]
    assert [p["compound_id"] for p in reaction["products"]] == ["m:b"]
    assert reaction["raw_reactants"][0]["compound_id"] == "m:b"


def test_map_snapshot_includes_isolated_catalog_records(tmp_path):
    source = tmp_path / "networkdb.json"
    source.write_text(json.dumps({
        "compounds": [{"id": "m:a"}, {"id": "m:b"}, {"id": "c:isolated"}],
        "reactions": [{"id": "r:1", "reactants": [{"compound_id": "m:a"}], "products": [{"compound_id": "m:b"}], "candidate_proteins": []}],
        "coverage": {"compound_records": 3},
    }))
    output = tmp_path / "map.json"
    result = build_map_snapshot(source, output)
    assert result["compounds"] == 3
    payload = json.loads(output.read_text())
    assert [c["id"] for c in payload["compounds"]] == ["m:a", "m:b", "c:isolated"]
    assert payload["focus"] == {"co2_reachable_compounds": 0, "reaction_connected_compounds": 2, "all_inventory_compounds": 3}
