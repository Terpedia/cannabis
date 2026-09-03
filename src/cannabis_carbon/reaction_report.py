from __future__ import annotations

import gzip
import json
from pathlib import Path

from .atom_mapping import map_reaction_smiles
from .terpedia import load_network


def build_reaction_report(source: Path, output: Path) -> dict:
    network = load_network(source)
    reactions = [e for e in network["entities"] if e.get("type") == "biochemical_reaction"]
    rows = []
    for reaction in reactions:
        attrs = reaction.get("attributes", {})
        mapping = map_reaction_smiles(attrs.get("reactionSmiles"))
        rows.append({"reaction_id": reaction["id"], "rhea_url": reaction.get("url"), "equation": attrs.get("equation"), "reaction_smiles": attrs.get("reactionSmiles"), **mapping})
    mapped = sum(sum(m["status"] == "inferred" for m in row["mappings"]) for row in rows)
    ambiguous = sum(sum(m["status"] == "ambiguous" for m in row["mappings"]) for row in rows)
    unresolved = sum(sum(m["status"] == "unresolved" for m in row["mappings"]) for row in rows)
    candidate = sum(sum(m["status"] == "candidate" for m in row["mappings"]) for row in rows)
    product_carbons = sum(row["product_carbon_atom_count"] for row in rows)
    report = {"schema": "cannabis-carbon.reaction-carbon-mapping.v1", "source": str(source), "reaction_count": len(rows), "status_counts": {status: sum(row["status"] == status for row in rows) for status in ("inferred", "unresolved")}, "carbon_counts": {"product_carbon_atoms": product_carbons, "mapped_carbon_atoms": mapped, "candidate_carbon_atoms": candidate, "mapped_carbon_atoms_including_candidates": mapped + candidate, "ambiguous_carbon_atoms": ambiguous, "unresolved_carbon_atoms": unresolved, "unresolved_or_ambiguous_carbon_atoms": ambiguous + unresolved, "mapping_coverage_percent": round(100 * mapped / product_carbons, 4) if product_carbons else 0.0}, "reactions": rows, "claim_boundary": "Structural atom mapping is not isotope tracing and does not establish in-vivo flux; CO2-transfer candidate assignments are not confirmed atom tracing."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {"reaction_count": len(rows), "status_counts": report["status_counts"]}
