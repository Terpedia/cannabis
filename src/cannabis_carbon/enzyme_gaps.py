"""Evidence-preserving audit of reactions lacking characterized enzymes."""

from __future__ import annotations

import json
from pathlib import Path


def build_enzyme_gap_audit(networkdb_path: Path, output: Path) -> dict:
    networkdb = json.loads(networkdb_path.read_text())
    gaps = []
    for reaction in networkdb.get("reactions", []):
        if reaction.get("enzyme_ids") or reaction.get("status") == "non_enzymatic":
            continue
        candidates = []
        for protein in reaction.get("candidate_proteins", []):
            sequence = protein.get("sequence_search") or {}
            diamond = protein.get("diamond_search") or {}
            candidates.append({"protein_id": protein.get("proteinId"), "accession": protein.get("accession"), "label": protein.get("label"), "candidate_origin": protein.get("candidateOrigin"), "exact_ec_numbers": protein.get("exactEcNumbers", []), "sequence_present": bool(sequence.get("sequence_present")), "sequence_length": sequence.get("length"), "sequence_sha256": sequence.get("sha256"), "diamond_hit_count": len(diamond.get("hits", [])), "strong_candidate_hit": bool(diamond.get("strong_candidate_hit")), "source_url": protein.get("sourceUrl")})
        gaps.append({"reaction_id": reaction.get("id"), "label": reaction.get("label"), "equation": reaction.get("equation"), "reaction_smiles": reaction.get("reaction_smiles"), "status": reaction.get("status"), "carbon_mapping_status": (reaction.get("carbon_mapping") or {}).get("status"), "source_url": reaction.get("source_url"), "candidate_protein_count": len(candidates), "candidates_with_sequence": sum(c["sequence_present"] for c in candidates), "candidates_with_diamond_hits": sum(c["diamond_hit_count"] > 0 for c in candidates), "strong_candidate_count": sum(c["strong_candidate_hit"] for c in candidates), "candidate_proteins": candidates, "claim_boundary": "These are genome-search and annotation candidates for an enzyme gap. They do not establish catalytic activity, substrate specificity, localization, or in-vivo Cannabis flux."})
    gaps.sort(key=lambda row: row["reaction_id"] or "")
    result = {"schema": "cannabis-carbon.enzyme-gap-audit.v1", "source_networkdb": str(networkdb_path), "reaction_gap_count": len(gaps), "reaction_gaps_with_candidates": sum(row["candidate_protein_count"] > 0 for row in gaps), "reaction_gaps_without_candidates": sum(row["candidate_protein_count"] == 0 for row in gaps), "candidate_protein_count": sum(row["candidate_protein_count"] for row in gaps), "candidate_proteins_with_sequence": sum(row["candidates_with_sequence"] for row in gaps), "candidate_proteins_with_diamond_hits": sum(row["candidates_with_diamond_hits"] for row in gaps), "strong_candidate_count": sum(row["strong_candidate_count"] for row in gaps), "reactions": gaps, "claim_boundary": "This audit identifies candidates for missing enzyme associations; only direct biochemical characterization can promote a candidate to a confirmed enzyme edge."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    return {key: result[key] for key in ("reaction_gap_count", "reaction_gaps_with_candidates", "candidate_protein_count", "strong_candidate_count")}
