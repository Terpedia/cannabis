"""Build reversible-upper-bound pathway hypotheses for candidate identity bridges."""

from __future__ import annotations

import json
from collections import defaultdict, deque
from pathlib import Path


def build_reversible_candidate_lineage(bridges_path: Path, lineage_path: Path, output: Path) -> dict:
    bridges = json.loads(bridges_path.read_text())
    lineage = json.loads(lineage_path.read_text())
    adjacency = defaultdict(list)
    for edge in lineage.get("carbon_edges", []):
        source, target = edge.get("reactant_entity_id"), edge.get("product_entity_id")
        if not source or not target:
            continue
        record = {"from": source, "to": target, "reaction_id": edge.get("reaction_id"), "status": edge.get("status"), "provenance": edge.get("provenance")}
        adjacency[source].append(record)
        adjacency[target].append({**record, "from": target, "to": source, "status": "reversible-upper-bound"})
    seed = lineage.get("co2_entity_id") or "chebi:16526"
    previous, previous_edge = {seed: None}, {}
    queue = deque([seed])
    while queue:
        node = queue.popleft()
        for edge in adjacency[node]:
            if edge["to"] in previous:
                continue
            previous[edge["to"]] = node
            previous_edge[edge["to"]] = edge
            queue.append(edge["to"])

    def path_to(node):
        path = []
        while node != seed and node in previous:
            path.append(previous_edge[node])
            node = previous[node]
        return list(reversed(path)) if node == seed else None

    rows = []
    for bridge in bridges.get("bridges", []):
        anchor = next((candidate for candidate in (bridge.get("core_precursor_entity_id"), bridge.get("core_product_entity_id")) if candidate in previous), None)
        path = path_to(anchor) if anchor else None
        if path is None:
            continue
        rows.append({"candidate_product_terpene_id": bridge.get("product_terpene_id"), "candidate_precursor_terpene_id": bridge.get("precursor_terpene_id"), "reaction_id": bridge.get("reaction_id"), "source_type": bridge.get("source_type"), "source_url": bridge.get("source_url"), "source_uniprot_id": bridge.get("source_uniprot_id"), "core_anchor_entity_id": anchor, "core_path_reaction_ids": [edge.get("reaction_id") for edge in path if edge.get("reaction_id")], "core_path_entity_ids": [seed] + [edge["to"] for edge in path], "core_path_step_count": len(path), "path_mode": "all-reactions-reversible-upper-bound", "status": "candidate", "claim_boundary": "This ordered path combines reversible structural reachability with a candidate identity bridge. It is a sensitivity hypothesis only; direction, exact identity, enzyme activity, isotope tracing, and in-vivo Cannabis production remain unestablished."})
    result = {"schema": "cannabis-carbon.terpene-identity-set-reversible-candidate-lineage.v1", "source_bridges": str(bridges_path), "source_lineage": str(lineage_path), "co2_entity_id": seed, "path_count": len(rows), "distinct_candidate_products": len({row["candidate_product_terpene_id"] for row in rows}), "distinct_core_anchors": len({row["core_anchor_entity_id"] for row in rows}), "rows": rows, "claim_boundary": "These are ordered reversible-upper-bound pathway hypotheses from CO2 to candidate identity-set structures. They are not directed biological pathways or proof of endogenous Cannabis biosynthesis."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    return {key: result[key] for key in ("path_count", "distinct_candidate_products", "distinct_core_anchors")}
