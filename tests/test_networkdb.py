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
