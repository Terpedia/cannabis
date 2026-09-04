import json
import gzip

from cannabis_carbon.inventory import build_specialty_inventory


def test_specialty_inventory_preserves_identity_and_reaction_status(tmp_path):
    compounds = tmp_path / "compounds.json"
    compounds.write_text(json.dumps({"compounds": [{"id": "CDB1", "label": "Cannabigerol", "aliases": [], "formula": "C21H32O2", "smiles": "CC", "carbon_atom_count": 21, "carbon_status": {"0": "unresolved"}, "source": "CDB", "source_url": "https://example.test"}]}))
    crosswalk = tmp_path / "crosswalk.json"
    crosswalk.write_text(json.dumps({"matches": [], "candidate_matches": []}))
    network = tmp_path / "network.json.gz"
    with gzip.open(network, "wt") as handle:
        json.dump({"entities": [], "statements": []}, handle)
    output = tmp_path / "inventory.json"
    result = build_specialty_inventory(compounds, crosswalk, network, output)
    assert result["record_count"] == 1
    assert json.loads(output.read_text())["records"][0]["identity_status"] == "unresolved"
