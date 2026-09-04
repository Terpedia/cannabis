from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path


def load_network(path: Path, additions_path: Path | None = None) -> dict:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        network = json.load(handle)
    addition_paths = [additions_path] if additions_path else [
        path.parent / "reaction-additions.json",
        path.parent / "varin-reaction-additions.json",
    ]
    for addition_path in addition_paths:
        if not addition_path.exists():
            continue
        additions = json.loads(addition_path.read_text())
        existing_entities = {entity["id"] for entity in network["entities"]}
        existing_statements = {(s.get("subjectId"), s.get("predicate"), s.get("objectEntityId")) for s in network["statements"]}
        for entity in additions.get("entities", []):
            if entity["id"] not in existing_entities:
                network["entities"].append(entity)
        for statement in additions.get("statements", []):
            key = (statement.get("subjectId"), statement.get("predicate"), statement.get("objectEntityId"))
            if key not in existing_statements:
                network["statements"].append(statement)
    return network


def cytoscape_elements(network: dict, direction_overrides: dict | None = None) -> dict:
    """Export compound-to-compound reaction edges with curated directions.

    Raw Rhea master records are often oriented opposite to the physiological
    direction.  The override layer is evidence, not a flux assertion; the raw
    participant orientation remains available in the edge metadata.
    """
    entities = {row["id"]: row for row in network["entities"]}
    reactions = {key: row for key, row in entities.items() if row.get("type") == "biochemical_reaction"}
    metabolites = {key: row for key, row in entities.items() if row.get("type") == "metabolite"}
    statements = network["statements"]
    reactants = defaultdict(list)
    products = defaultdict(list)
    enzymes = defaultdict(list)
    evidence = defaultdict(list)
    for statement in statements:
        predicate = statement["predicate"]
        subject = statement["subjectId"]
        obj = statement["objectEntityId"]
        if predicate == "has_reactant" and subject in reactions:
            reactants[subject].append(obj)
        elif predicate == "has_product" and subject in reactions:
            products[subject].append(obj)
        elif predicate in ("catalyzes", "maps_to_reaction", "has_catalytic_activity") and obj in reactions:
            enzymes[obj].append(subject)
            evidence[obj].append(statement.get("qualifiers", {}))
    nodes = []
    for metabolite_id, row in metabolites.items():
        nodes.append({"data": {"id": metabolite_id, "label": row.get("label", metabolite_id), "kind": "compound", "status": "supported", "source_url": row.get("url")}})
    edges = []
    for reaction_id in reactions:
        direction = (direction_overrides or {}).get(reaction_id, {})
        oriented_reactants = products[reaction_id] if direction.get("orientation") == "reverse_master" else reactants[reaction_id]
        oriented_products = reactants[reaction_id] if direction.get("orientation") == "reverse_master" else products[reaction_id]
        enzyme_ids = sorted(set(enzymes[reaction_id]))
        enzyme_labels = [entities[e].get("label", e) for e in enzyme_ids if e in entities]
        reaction_class = reactions[reaction_id].get("attributes", {}).get("reactionClass", "")
        status = "non_enzymatic" if str(reaction_class).startswith("non-enzymatic-") else "supported" if any(q.get("directExperimentalEvidence") for q in evidence[reaction_id]) else "candidate"
        label = reactions[reaction_id].get("label", reaction_id)
        if enzyme_labels:
            label = f"{label} · {'; '.join(enzyme_labels[:3])}"
        else:
            label = f"{label} · enzyme unresolved"
        for reactant in oriented_reactants:
            for product in oriented_products:
                if reactant in metabolites and product in metabolites:
                    edges.append({"data": {
                        "id": f"{reaction_id}:{reactant}>{product}",
                        "source": reactant,
                        "target": product,
                        "label": label,
                        "kind": "reaction",
                        "status": status,
                        "reaction_id": reaction_id,
                        "enzyme_ids": enzyme_ids,
                        "source_url": reactions[reaction_id].get("url"),
                        "raw_direction": "product-to-reactant" if direction.get("orientation") == "reverse_master" else "reactant-to-product",
                        "directional_rhea_id": direction.get("directional_rhea_id"),
                        "direction_evidence_url": direction.get("source"),
                        "direction_reason": direction.get("reason"),
                    }})
    return {"schema": "cannabis-carbon.cytoscape.v1", "nodes": nodes, "edges": edges, "stats": {"metabolites": len(nodes), "reaction_edges": len(edges), "reactions": len(reactions), "direction_overrides": sum(bool((direction_overrides or {}).get(rid)) for rid in reactions)}}
