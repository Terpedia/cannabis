import json

from cannabis_carbon.hypotheses import build_candidate_queue


def test_candidate_queue_preserves_missing_producer_boundary(tmp_path):
    source = tmp_path / "hypotheses.json"
    source.write_text(json.dumps({"summary": {"coverageBoundary": "not proof"}, "hypotheses": {"missingProducers": [{"hypothesisId": "h1", "hypothesisType": "missing-producing-reaction", "status": "structurally_possible", "targetMetabolite": {"entityId": "chebi:1"}, "candidateReaction": {}, "candidateEnzymes": [], "claimBoundary": "unknown"}]}}))
    output = tmp_path / "queue.json"
    summary = build_candidate_queue(source, output)
    assert summary["missing_producer"] == 1
    assert json.loads(output.read_text())["items"][0]["status"] == "structurally_possible"


def test_candidate_queue_merges_specialized_enzyme_search(tmp_path):
    source = tmp_path / "hypotheses.json"
    source.write_text(json.dumps({"summary": {"coverageBoundary": "not proof"}, "hypotheses": {"missingProducers": [{"hypothesisId": "h1", "hypothesisType": "missing-producing-reaction", "status": "candidate", "targetMetabolite": {}, "candidateReaction": {"reactionId": "rhea:1"}, "candidateEnzymes": [], "claimBoundary": "unknown"}]}}))
    (tmp_path / "enzyme-candidate-additions.json").write_text(json.dumps({"searches": [{"reaction_id": "rhea:1", "candidate_proteins": [{"proteinId": "uniprot:p1", "accession": "P1"}]}]}))
    output = tmp_path / "queue.json"
    build_candidate_queue(source, output)
    item = json.loads(output.read_text())["items"][0]
    assert item["candidate_proteins"][0]["proteinId"] == "uniprot:p1"


def test_target_hypothesis_does_not_assert_catalog_presence_is_endogenous(tmp_path):
    from cannabis_carbon.test_hypotheses import _target_hypotheses

    lineage = tmp_path / "lineage.json"
    lineage.write_text(json.dumps({"targets": [{"cannabisdb_id": "CDB1", "status": "unresolved", "reason": "no-path"}]}))
    compounds = tmp_path / "compounds.json"
    compounds.write_text(json.dumps({"compounds": [{"id": "CDB1", "label": "test", "aliases": [], "carbon_atom_count": 1}]}))
    target = _target_hypotheses(lineage, compounds)[0]
    assert target["claim"].startswith("Test whether CDB1")
