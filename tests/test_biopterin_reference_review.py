import hashlib
import json
from pathlib import Path


def test_generic_biopterin_lead_preserves_annotation_and_gap_boundary():
    root = Path(__file__).resolve().parents[1]
    review = json.loads((root / "data/curation/biopterin-reference-review.json").read_text())
    raw = (root / review["source"]["raw_path"]).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == review["source"]["sha256"]
    protein = json.loads(raw)
    assert protein["primaryAccession"] == review["source"]["accession"]
    activities = [c for c in protein["comments"] if c["commentType"] == "CATALYTIC ACTIVITY"]
    for observation in review["annotation_observations"]:
        activity = next(c for c in activities if any(
            x["id"] == observation["rhea_id"] for x in c["reaction"]["reactionCrossReferences"]
        ))
        assert activity["reaction"]["ecNumber"] == observation["ec"]
        assert {f'{e["source"]}:{e["id"]}' for e in activity["reaction"]["evidences"]} == set(observation["activity_evidence"])
        assert all(e["evidenceCode"] == observation["activity_evidence_code"] for e in activity["reaction"]["evidences"])
        direction = activity["physiologicalReactions"][0]
        assert direction["directionType"] == observation["direction_in_source_equation"]
        assert direction["reactionCrossReference"]["id"] == observation["physiological_rhea_id"]
        assert all(e["evidenceCode"] == observation["direction_evidence_code"] for e in direction["evidences"])
        assert observation["rhea_id"] not in review["exact_rhea_family"]
    parent = json.loads((root / "data/reports/phase1-nonplant-reference-review.json").read_text())
    row = next(r for r in parent["rows"] if r["reaction_id"] == review["reaction_id"])
    assert row["target_ids"] == review["target_ids"]
    assert set(next(iter(row["rhea_families"].values())).values()) == set(review["exact_rhea_family"])
    assert review["interpretation"]["exact_reaction_annotation_match"] is False
    assert review["interpretation"]["candidate_model_changed"] is False
    assert review["interpretation"]["genome_screen_performed"] is False
