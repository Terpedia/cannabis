"""Candidate CO2 reachability through source-linked hypothesis edges."""

from __future__ import annotations

import json
from collections import Counter, deque
from pathlib import Path


def build_hypothesis_lineage(networkdb_path: Path, output: Path) -> dict:
    """Traverse candidate hypothesis edges without promoting them to core lineage.

    Core CO2 reachability is the seed.  Only explicitly ``candidate`` edges are
    traversed; unresolved-substrate edges remain visible as blockers but never
    create reachability.  This produces a falsifiable candidate layer for
    target triage rather than a claim of confirmed biosynthesis.
    """
    networkdb = json.loads(networkdb_path.read_text())
    compounds = networkdb.get("compounds", [])
    edges = networkdb.get("hypothetical_connections", [])
    by_id = {compound.get("id"): compound for compound in compounds}
    seeds = {
        compound["id"]
        for compound in compounds
        if (compound.get("co2_reachable_carbon_atoms") or 0) > 0
    }
    candidate_edges = [edge for edge in edges if edge.get("status") == "candidate"]
    outgoing: dict[str, list[dict]] = {}
    for edge in candidate_edges:
        outgoing.setdefault(edge.get("substrate_compound_id"), []).append(edge)

    reachable = set(seeds)
    parent: dict[str, tuple[str, dict]] = {}
    queue = deque(seeds)
    while queue:
        source = queue.popleft()
        for edge in outgoing.get(source, []):
            target = edge.get("product_compound_id")
            if target and target not in reachable:
                reachable.add(target)
                parent[target] = (source, edge)
                queue.append(target)

    target_rows = []
    for compound in compounds:
        if compound.get("namespace") != "cannabisdb":
            continue
        compound_id = compound["id"]
        core_atoms = compound.get("co2_reachable_carbon_atoms") or 0
        if compound_id not in reachable:
            status = "unresolved"
            reason = "no-core-or-candidate-hypothesis-path"
            path = []
        elif compound_id in seeds:
            status = "core"
            reason = "core-carbon-lineage"
            path = []
        else:
            status = "candidate"
            reason = "candidate-hypothesis-path-from-core-co2-lineage"
            path = []
            current = compound_id
            while current in parent:
                source, edge = parent[current]
                path.append({
                    "reaction_id": edge.get("reaction_id"),
                    "from_compound_id": source,
                    "to_compound_id": current,
                    "source_url": edge.get("source_url"),
                    "evidence_type": edge.get("evidence_type"),
                    "enzyme_evidence_count": len(edge.get("enzyme_evidence") or []),
                    "enzyme_catalog_count": len(edge.get("enzyme_catalog") or []),
                    "balance_status": edge.get("balance_status"),
                })
                current = source
            path.reverse()
        target_rows.append({
            "cannabisdb_id": compound_id,
            "label": compound.get("label"),
            "carbon_atom_count": compound.get("carbon_atom_count", 0),
            "core_reachable_carbon_atoms": core_atoms,
            "status": status,
            "reason": reason,
            "path": path,
        })

    status_counts = Counter(row["status"] for row in target_rows)
    carbon_counts = Counter()
    for row in target_rows:
        carbon_counts[row["status"]] += row["carbon_atom_count"]
    blocked_edges = Counter(edge.get("blocker") or "unspecified" for edge in edges if edge.get("status") == "unresolved")
    report = {
        "schema": "cannabis-carbon.hypothesis-lineage.v1",
        "source_networkdb": str(networkdb_path),
        "carbon_source_policy": "CO2 is the only admissible carbon source; candidate hypothesis paths are provisional and separate from the core balanced lineage.",
        "seed_core_reachable_compounds": len(seeds),
        "candidate_edges_traversed": len(candidate_edges),
        "candidate_reachable_entities": len(reachable - seeds),
        "target_summary": {
            "counts_by_status": dict(sorted(status_counts.items())),
            "carbon_atoms_by_status": dict(sorted(carbon_counts.items())),
            "candidate_cannabisdb_targets": sum(row["status"] == "candidate" for row in target_rows),
            "candidate_cannabisdb_carbon_atoms": carbon_counts["candidate"],
        },
        "blocked_unresolved_hypothesis_edges": dict(sorted(blocked_edges.items())),
        "targets": target_rows,
        "claim_boundary": "Candidate paths are source-linked structural hypotheses, not confirmed enzyme activity, isotope tracing, or proof of in-vivo Cannabis biosynthesis. Unresolved edges are not traversed.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return report
