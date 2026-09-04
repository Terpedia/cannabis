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
        json.dump({"entities": [], "statements": [{"subjectId": "rhea:1", "predicate": "has_product", "objectEntityId": "cannabisdb:CDB1", "qualifiers": {"stoichiometricCoefficient": 1}}]}, handle)
    output = tmp_path / "inventory.json"
    result = build_specialty_inventory(compounds, crosswalk, network, output)
    assert result["record_count"] == 1
    assert json.loads(output.read_text())["records"][0]["identity_status"] == "unresolved"
    assert result["records_with_reaction_participation"] == 1


def test_specialty_inventory_adds_carbon_weighted_review_queue(tmp_path):
    compounds = tmp_path / "compounds.json"
    compounds.write_text(json.dumps({"compounds": [{"id": "CDB1", "label": "Cannabigerol", "aliases": [], "formula": "C21H32O2", "smiles": "CC", "carbon_atom_count": 21}]}))
    crosswalk = tmp_path / "crosswalk.json"
    crosswalk.write_text(json.dumps({"matches": [], "candidate_matches": []}))
    network = tmp_path / "network.json.gz"
    with gzip.open(network, "wt") as handle:
        json.dump({"entities": [], "statements": []}, handle)
    lineage = tmp_path / "lineage.json"
    lineage.write_text(json.dumps({"targets": [{"cannabisdb_id": "CDB1", "status": "candidate", "reachable_carbon_atoms": 4, "reversible_upper_bound_reachable_carbon_atoms": 9, "reason": "external-carbon-input"}]}))
    pubchem = tmp_path / "pubchem.json"
    pubchem.write_text(json.dumps({"records": [{"cannabisdb_id": "CDB1", "status": "resolved", "pubchem": {"CID": 123}}]}))
    output = tmp_path / "inventory.json"
    build_specialty_inventory(compounds, crosswalk, network, output, lineage, pubchem)
    report = json.loads(output.read_text())
    record = report["records"][0]
    assert record["unresolved_carbon_atoms"] == 17
    assert record["reversible_upper_bound_reachable_carbon_atoms"] == 9
    assert record["pubchem_cid"] == 123
    assert report["review_queue"][0]["unresolved_carbon_atoms"] == 17
