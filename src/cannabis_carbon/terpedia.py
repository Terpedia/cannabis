from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path


def load_network(path: Path, additions_path: Path | None = None) -> dict:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        network = json.load(handle)
    additions_path = additions_path or path.parent / "reaction-additions.json"
    if additions_path.exists():
        additions = json.loads(additions_path.read_text())
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


def cytoscape_elements(network: dict) -> dict:
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
        enzyme_ids = sorted(set(enzymes[reaction_id]))
        enzyme_labels = [entities[e].get("label", e) for e in enzyme_ids if e in entities]
        status = "supported" if any(q.get("directExperimentalEvidence") for q in evidence[reaction_id]) else "candidate"
        label = reactions[reaction_id].get("label", reaction_id)
        if enzyme_labels:
            label = f"{label} · {'; '.join(enzyme_labels[:3])}"
        else:
            label = f"{label} · enzyme unresolved"
        for reactant in reactants[reaction_id]:
            for product in products[reaction_id]:
                if reactant in metabolites and product in metabolites:
                    edges.append({"data": {"id": f"{reaction_id}:{reactant}>{product}", "source": reactant, "target": product, "label": label, "kind": "reaction", "status": status, "reaction_id": reaction_id, "enzyme_ids": enzyme_ids, "source_url": reactions[reaction_id].get("url")}})
    return {"schema": "cannabis-carbon.cytoscape.v1", "nodes": nodes, "edges": edges, "stats": {"metabolites": len(nodes), "reaction_edges": len(edges), "reactions": len(reactions)}}
