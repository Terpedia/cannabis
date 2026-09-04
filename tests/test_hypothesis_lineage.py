import json

from cannabis_carbon.hypothesis_lineage import build_hypothesis_lineage


def test_candidate_lineage_traverses_only_candidate_edges(tmp_path):
    networkdb = tmp_path / "networkdb.json"
    networkdb.write_text(json.dumps({
        "compounds": [
            {"id": "chebi:16526", "namespace": "terpedia", "co2_reachable_carbon_atoms": 1, "carbon_atom_count": 1},
            {"id": "CDB1", "namespace": "cannabisdb", "co2_reachable_carbon_atoms": 0, "carbon_atom_count": 2},
            {"id": "CDB2", "namespace": "cannabisdb", "co2_reachable_carbon_atoms": 0, "carbon_atom_count": 3},
        ],
        "hypothetical_connections": [
            {"substrate_compound_id": "chebi:16526", "product_compound_id": "CDB1", "status": "candidate", "reaction_id": "R1"},
            {"substrate_compound_id": "CDB1", "product_compound_id": "CDB2", "status": "unresolved", "blocker": "missing-corpus-substrate", "reaction_id": "R2"},
        ],
    }))
    report = build_hypothesis_lineage(networkdb, tmp_path / "lineage.json")
    assert report["target_summary"]["counts_by_status"] == {"candidate": 1, "unresolved": 1}
    assert report["blocked_unresolved_hypothesis_edges"] == {"missing-corpus-substrate": 1}
