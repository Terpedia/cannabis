from __future__ import annotations

import gzip
import json
from pathlib import Path


def load_hypotheses(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_candidate_queue(path: Path, output: Path) -> dict:
    """Normalize Terpedia's unresolved pathway hypotheses for genome searches."""
    source = load_hypotheses(path)
    queue = []
    for row in source["hypotheses"].get("unresolvedExactEc", []):
        queue.append({"id": row["hypothesisId"], "kind": row["hypothesisType"], "status": row["status"], "reaction_id": row.get("reactionId"), "ec_number": row.get("ecNumber"), "target_metabolite": None, "candidate_proteins": row.get("candidateEnzymes", []), "source": row.get("claimBoundary")})
    for row in source["hypotheses"].get("unresolvedPartialEc", []):
        queue.append({"id": row["proteinId"], "kind": row["hypothesisReason"], "status": "candidate", "reaction_id": None, "ec_number": (row.get("partialEcNumbers") or [None])[0], "target_metabolite": None, "candidate_proteins": [row], "source": row.get("sourceUrl")})
    for row in source["hypotheses"].get("missingProducers", []):
        target = row.get("targetMetabolite", {})
        queue.append({"id": row["hypothesisId"], "kind": row["hypothesisType"], "status": row["status"], "reaction_id": row.get("candidateReaction", {}).get("reactionId"), "ec_number": None, "target_metabolite": target, "candidate_proteins": row.get("candidateEnzymes", []), "downstream_reactions": row.get("downstreamConsumers", []), "source": row.get("claimBoundary")})
    for row in source["hypotheses"].get("blockedKnownReactions", []):
        reaction = row.get("reaction", {})
        target = row.get("targetMetabolite", {})
        queue.append({"id": row["hypothesisId"], "kind": row["hypothesisType"], "status": row["status"], "reaction_id": reaction.get("reactionId"), "ec_number": None, "target_metabolite": target, "candidate_proteins": row.get("candidateEnzymes", []), "source": row.get("claimBoundary")})
    addition_paths = [path.parent / "reaction-additions.json", path.parent / "varin-reaction-additions.json"]
    for additions_path in addition_paths:
        if not additions_path.exists():
            continue
        additions = json.loads(additions_path.read_text())
        known_reactions = {item.get("reaction_id") for item in queue}
        for reaction in additions.get("entities", []):
            if reaction.get("type") != "biochemical_reaction" or reaction.get("id") in known_reactions:
                continue
            queue.append({"id": f"missing-producer:{reaction['id']}", "kind": "curated-reaction-addition-missing-producer", "status": "candidate", "reaction_id": reaction["id"], "ec_number": (reaction.get("attributes", {}).get("ecNumbers") or [None])[0], "target_metabolite": None, "candidate_proteins": [], "source": reaction.get("source_url") or reaction.get("url")})
    candidate_additions_path = path.parent / "enzyme-candidate-additions.json"
    if candidate_additions_path.exists():
        candidate_additions = json.loads(candidate_additions_path.read_text())
        items_by_reaction = {item.get("reaction_id"): item for item in queue}
        for search in candidate_additions.get("searches", []):
            item = items_by_reaction.get(search.get("reaction_id"))
            if item is None:
                continue
            existing = {protein.get("proteinId") for protein in item["candidate_proteins"]}
            for protein in search.get("candidate_proteins", []):
                if protein.get("proteinId") not in existing:
                    item["candidate_proteins"].append(protein)
    report = {"schema": "cannabis-carbon.candidate-work-queue.v1", "source": str(path), "summary": {"total": len(queue), "with_candidate_proteins": sum(bool(row["candidate_proteins"]) for row in queue), "missing_producer": len(source["hypotheses"].get("missingProducers", [])), "unresolved_ec": len(source["hypotheses"].get("unresolvedExactEc", [])) + len(source["hypotheses"].get("unresolvedPartialEc", [])), "blocked_reaction": len(source["hypotheses"].get("blockedKnownReactions", []))}, "items": queue, "claim_boundary": source["summary"]["coverageBoundary"]}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return report["summary"]
