import gzip
import json

from cannabis_carbon.networkdb import build_networkdb


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
