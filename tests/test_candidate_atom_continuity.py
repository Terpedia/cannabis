import json

from cannabis_carbon.candidate_lineage import atom_continuity_blockers, build_reversible_candidate_lineage


def test_connected_entities_do_not_imply_connected_atoms():
    edges = [
        {"from_entity_id": "co2", "from_atom": 1, "to_entity_id": "a", "to_atom": 3},
        {"from_entity_id": "a", "from_atom": 7, "to_entity_id": "b", "to_atom": 0},
    ]
    assert atom_continuity_blockers(edges)[0]["reason"] == "carbon-atom-discontinuity"
    edges[1]["from_atom"] = 3
    assert atom_continuity_blockers(edges) == []


def test_reverse_traversal_swaps_atom_indices(tmp_path):
    bridges = tmp_path / "bridges.json"
    lineage = tmp_path / "lineage.json"
    output = tmp_path / "out.json"
    bridges.write_text(json.dumps({"bridges": [{"core_precursor_entity_id": "a"}]}))
    lineage.write_text(json.dumps({"co2_entity_id": "co2", "carbon_edges": [{
        "reactant_entity_id": "a", "reactant_atom": 15,
        "product_entity_id": "co2", "product_atom": 1,
    }]}))
    build_reversible_candidate_lineage(bridges, lineage, output)
    row = json.loads(output.read_text())["rows"][0]
    assert row["core_path_carbon_edges"][0]["from_atom"] == 1
    assert row["core_path_carbon_edges"][0]["to_atom"] == 15
    assert row["core_atom_continuity"]["status"] == "continuous"
    assert row["carbon_provenance_status"] == "unresolved"
