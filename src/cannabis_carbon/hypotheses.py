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
    additions_path = path.parent / "reaction-additions.json"
    if additions_path.exists():
        additions = json.loads(additions_path.read_text())
        known_reactions = {item.get("reaction_id") for item in queue}
        for reaction in additions.get("entities", []):
            if reaction.get("type") != "biochemical_reaction" or reaction.get("id") in known_reactions:
                continue
            queue.append({"id": f"missing-producer:{reaction['id']}", "kind": "curated-reaction-addition-missing-producer", "status": "candidate", "reaction_id": reaction["id"], "ec_number": (reaction.get("attributes", {}).get("ecNumbers") or [None])[0], "target_metabolite": None, "candidate_proteins": [], "source": reaction.get("source_url") or reaction.get("url")})
    report = {"schema": "cannabis-carbon.candidate-work-queue.v1", "source": str(path), "summary": {"total": len(queue), "with_candidate_proteins": sum(bool(row["candidate_proteins"]) for row in queue), "missing_producer": len(source["hypotheses"].get("missingProducers", [])), "unresolved_ec": len(source["hypotheses"].get("unresolvedExactEc", [])) + len(source["hypotheses"].get("unresolvedPartialEc", [])), "blocked_reaction": len(source["hypotheses"].get("blockedKnownReactions", []))}, "items": queue, "claim_boundary": source["summary"]["coverageBoundary"]}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return report["summary"]
