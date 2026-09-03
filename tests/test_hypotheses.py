import json

from cannabis_carbon.hypotheses import build_candidate_queue


def test_candidate_queue_preserves_missing_producer_boundary(tmp_path):
    source = tmp_path / "hypotheses.json"
    source.write_text(json.dumps({"summary": {"coverageBoundary": "not proof"}, "hypotheses": {"missingProducers": [{"hypothesisId": "h1", "hypothesisType": "missing-producing-reaction", "status": "structurally_possible", "targetMetabolite": {"entityId": "chebi:1"}, "candidateReaction": {}, "candidateEnzymes": [], "claimBoundary": "unknown"}]}}))
    output = tmp_path / "queue.json"
    summary = build_candidate_queue(source, output)
    assert summary["missing_producer"] == 1
    assert json.loads(output.read_text())["items"][0]["status"] == "structurally_possible"
