import json

from cannabis_carbon.candidate_balance import audit_candidate_expansion_balances


def test_candidate_expansion_balance_deduplicates_reaction_edges(tmp_path):
    source = tmp_path / "expansion.json"
    source.write_text(json.dumps({"rows": [
        {"reaction_id": "R1", "reaction_smarts": "CC>>CC", "source_type": "Rhea"},
        {"reaction_id": "R1", "reaction_smarts": "CC>>CC", "source_type": "MARTS-DB"},
        {"reaction_id": "R2", "reaction_smarts": "CC>>C", "source_type": "Rhea"},
    ]}))
    output = tmp_path / "audit.json"
    summary = audit_candidate_expansion_balances(source, output)
    assert summary["edge_count"] == 3
    assert summary["unique_reaction_count"] == 2
    assert summary["edge_status_counts"] == {"balanced": 2, "imbalanced": 1}
    report = json.loads(output.read_text())
    assert report["reactions"][0]["edge_count"] == 2
    assert report["reactions"][0]["source_types"] == ["MARTS-DB", "Rhea"]
