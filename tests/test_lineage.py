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
