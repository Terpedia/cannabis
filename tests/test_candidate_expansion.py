import gzip
import json

from cannabis_carbon.identity_set_paths import build_candidate_expansion_bridges, build_candidate_expansion_carbon_mapping


def test_candidate_expansion_bridge_preserves_reachability_modes(tmp_path):
    expansion = tmp_path / "expansion.json"
    expansion.write_text(json.dumps({"rows": [{
        "product_terpene_id": "T2", "precursor_terpene_id": "T1",
        "product_smiles": "CCC", "precursor_smiles": "CC",
        "reaction_id": "R1", "source_type": "Rhea", "expansion_depth": 1,
    }]}))
    network = tmp_path / "network.json.gz"
    with gzip.open(network, "wt") as handle:
        json.dump({"entities": [
            {"id": "m:core-product", "type": "metabolite", "attributes": {"canonicalSmiles": "CCC"}},
            {"id": "m:core-precursor", "type": "metabolite", "attributes": {"canonicalSmiles": "CC"}},
        ], "statements": []}, handle)
    lineage = tmp_path / "lineage.json"
    lineage.write_text(json.dumps({"reachable_carbon_entity_ids": [], "reversible_upper_bound_reachable_carbon_entity_ids": ["m:core-precursor"]}))
    output = tmp_path / "bridges.json"
    result = build_candidate_expansion_bridges(expansion, network, lineage, output)
    payload = json.loads(output.read_text())
    assert result["bridge_count"] == 1
    assert payload["bridges"][0]["touches_co2_reachable_core"] is False
    assert payload["bridges"][0]["touches_reversible_co2_core"] is True


def test_candidate_expansion_carbon_mapping_retains_unresolved_carbons(tmp_path):
    expansion = tmp_path / "expansion.json"
    expansion.write_text(json.dumps({"rows": [{
        "product_terpene_id": "T2", "precursor_terpene_id": "T1",
        "product_smiles": "CCCC", "precursor_smiles": "CC",
        "reaction_id": "R1", "source_type": "Rhea", "expansion_depth": 1,
    }]}))
    bridges = tmp_path / "bridges.json"
    bridges.write_text(json.dumps({"bridges": [{
        "product_terpene_id": "T2", "precursor_terpene_id": "T1",
        "reaction_id": "R1", "expansion_depth": 1,
    }]}))
    output = tmp_path / "mapping.json"
    result = build_candidate_expansion_carbon_mapping(expansion, bridges, output)
    assert result["pair_count"] == 1
    assert result["unresolved_product_carbon_atoms"] > 0
    assert json.loads(output.read_text())["rows"][0]["status"] == "unresolved"
