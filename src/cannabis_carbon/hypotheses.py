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


def build_carbon_mapping_queue(mapping_path: Path, networkdb_path: Path, output: Path) -> dict:
    """Rank reaction carbon-mapping blockers for atom-lineage curation."""
    mapping = json.loads(mapping_path.read_text())
    networkdb = json.loads(networkdb_path.read_text())
    reaction_by_id = {reaction["id"]: reaction for reaction in networkdb.get("reactions", [])}
    items = []
    for row in mapping.get("reactions", []):
        counts = {status: sum(mapping_row.get("status") == status for mapping_row in row.get("mappings", [])) for status in ("inferred", "candidate", "ambiguous", "unresolved")}
        blocked = counts["ambiguous"] + counts["unresolved"]
        if not blocked:
            continue
        reaction = reaction_by_id.get(row["reaction_id"], {})
        items.append({
            "rank_key": (counts["unresolved"], counts["ambiguous"], blocked),
            "reaction_id": row["reaction_id"],
            "label": reaction.get("label", row["reaction_id"]),
            "equation": reaction.get("equation"),
            "source_url": reaction.get("source_url") or row.get("rhea_url"),
            "product_carbon_atom_count": row.get("product_carbon_atom_count", 0),
            "mapping_status_counts": counts,
            "blocker": "unresolved-or-ambiguous-rdkit-carbon-mapping",
            "enzyme_ids": reaction.get("enzyme_ids", []),
            "candidate_protein_ids": sorted({protein.get("proteinId") for protein in reaction.get("candidate_proteins", []) if protein.get("proteinId")}),
            "direction": reaction.get("direction", {}),
            "direction_status": "curated" if reaction.get("direction", {}).get("directional_rhea_id") else "raw",
            "mapping_methods": sorted({mapping_row.get("method") for mapping_row in row.get("mappings", []) if mapping_row.get("method")}),
        })
    items.sort(key=lambda item: tuple(-value for value in item["rank_key"]) + (item["reaction_id"],))
    for rank, item in enumerate(items, 1):
        item["rank"] = rank
        item.pop("rank_key")
    report = {
        "schema": "cannabis-carbon.carbon-mapping-work-queue.v1",
        "source_mapping": str(mapping_path),
        "source_networkdb": str(networkdb_path),
        "summary": {
            "reactions_with_mapping_blockers": len(items),
            "ambiguous_product_carbon_rows": sum(item["mapping_status_counts"]["ambiguous"] for item in items),
            "unresolved_product_carbon_rows": sum(item["mapping_status_counts"]["unresolved"] for item in items),
            "total_blocked_product_carbon_rows": sum(item["rank_key"] if "rank_key" in item else item["mapping_status_counts"]["ambiguous"] + item["mapping_status_counts"]["unresolved"] for item in items),
        },
        "items": items,
        "claim_boundary": "This queue prioritizes structural atom-mapping work. Resolving a mapping does not establish reaction direction, enzyme function, or in-vivo cannabis biosynthesis.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return report["summary"]
