"""Directed carbon-lineage audit from CO2 through the Terpedia reaction graph."""

from __future__ import annotations

import json
import re
from collections import defaultdict, deque
from functools import lru_cache
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger

from .atom_mapping import _molecules, map_reaction_smiles
from .terpedia import load_network

RDLogger.DisableLog("rdApp.*")


@lru_cache(maxsize=None)
def _key(smiles: str | None) -> str | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        return None
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return None


def _participant_ids(network: dict, reaction_id: str, predicate: str) -> list[str]:
    return [s["objectEntityId"] for s in network["statements"] if s.get("subjectId") == reaction_id and s.get("predicate") == predicate]


def _direction_evidence(network: dict) -> tuple[dict, list[str]]:
    directions = {}
    conflicts = []
    by_reaction = defaultdict(set)
    statement_sources = {}
    for statement in network["statements"]:
        predicate = statement.get("predicate")
        if predicate not in ("physiological_direction_left_to_right", "physiological_direction_right_to_left"):
            continue
        reaction_id = statement.get("subjectId")
        by_reaction[reaction_id].add(predicate)
        statement_sources[reaction_id] = (statement.get("sources") or [{}])[0].get("url")
    for reaction_id, predicates in by_reaction.items():
        if len(predicates) > 1:
            conflicts.append(reaction_id)
        elif "physiological_direction_right_to_left" in predicates:
            directions[reaction_id] = {"directional_rhea_id": None, "orientation": "reverse_master", "source": statement_sources.get(reaction_id), "reason": "Terpedia physiological direction statement"}
    return directions, conflicts


def _resolve_molecules(molecules, participant_ids, entities):
    by_key = defaultdict(list)
    for entity_id in participant_ids:
        key = _key(entities.get(entity_id, {}).get("attributes", {}).get("canonicalSmiles"))
        if key:
            by_key[key].append(entity_id)
    resolved = []
    for molecule in molecules:
        key = _key(Chem.MolToSmiles(molecule))
        candidates = by_key.get(key, [])
        resolved.append(candidates[0] if len(candidates) == 1 else None)
    return resolved


def _carbon_count(smiles: str | None) -> int:
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()) if mol is not None else 0


def _carbon_atom_indices(smiles: str | None) -> list[int]:
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6] if mol is not None else []


def build_carbon_lineage(network_path: Path, mapping_path: Path, crosswalk_path: Path, compounds_path: Path, output: Path, directions_path: Path | None = None) -> dict:
    network = load_network(network_path)
    mapping_report = json.loads(mapping_path.read_text())
    crosswalk = json.loads(crosswalk_path.read_text())
    compounds = json.loads(compounds_path.read_text())["compounds"]
    directions, direction_conflicts = _direction_evidence(network)
    if directions_path and directions_path.exists():
        directions.update(json.loads(directions_path.read_text()))
    entities = {e["id"]: e for e in network["entities"]}
    edges = []
    edge_blocks = defaultdict(int)
    for row in mapping_report["reactions"]:
        if not row.get("reaction_smiles"):
            edge_blocks["missing-reaction-smiles"] += 1
            continue
        direction = directions.get(row["reaction_id"], {})
        reaction_smiles = row["reaction_smiles"]
        reactant_predicate, product_predicate = "has_reactant", "has_product"
        if direction.get("orientation") == "reverse_master":
            left, right = reaction_smiles.split(">>", 1)
            reaction_smiles = f"{right}>>{left}"
            reactant_predicate, product_predicate = "has_product", "has_reactant"
        left, right = reaction_smiles.split(">>", 1)
        reactants, products = _molecules(left), _molecules(right)
        reaction_id = row["reaction_id"]
        reactant_ids = _resolve_molecules(reactants, _participant_ids(network, reaction_id, reactant_predicate), entities)
        product_ids = _resolve_molecules(products, _participant_ids(network, reaction_id, product_predicate), entities)
        mappings = map_reaction_smiles(reaction_smiles)["mappings"] if direction.get("orientation") == "reverse_master" else row["mappings"]
        for mapping in mappings:
            if mapping.get("status") not in ("inferred", "candidate", "ambiguous"):
                continue
            alternatives = mapping.get("alternatives") or [{"reactant_index": mapping.get("reactant_index"), "reactant_atom": mapping.get("reactant_atom")}]
            for alternative in alternatives:
                ri, pi = alternative["reactant_index"], mapping["product_index"]
                if ri is None or ri >= len(reactant_ids) or pi >= len(product_ids) or not reactant_ids[ri] or not product_ids[pi]:
                    edge_blocks["participant-structure-unresolved"] += 1
                    continue
                edges.append({"reaction_id": reaction_id, "reactant_entity_id": reactant_ids[ri], "reactant_atom": alternative["reactant_atom"], "product_entity_id": product_ids[pi], "product_atom": mapping["product_atom"], "status": "candidate" if mapping["status"] == "ambiguous" else mapping["status"], "provenance": direction.get("source") or row.get("rhea_url"), "directional_rhea_id": direction.get("directional_rhea_id"), "mapping_reason": mapping.get("reason")})

    forward = defaultdict(set)
    for edge in edges:
        forward[(edge["reactant_entity_id"], edge["reactant_atom"])].add(((edge["product_entity_id"], edge["product_atom"]), edge["status"]))
    co2_id = "chebi:16526"
    co2_smiles = entities.get(co2_id, {}).get("attributes", {}).get("canonicalSmiles")
    co2_atoms = [(co2_id, i) for i in _carbon_atom_indices(co2_smiles)]
    reachable = set(co2_atoms)
    candidate_reachable = set()
    seen_states = {(node, False) for node in co2_atoms}
    queue = deque((node, False) for node in co2_atoms)
    while queue:
        node, has_candidate = queue.popleft()
        for child, edge_status in forward[node]:
            child_has_candidate = has_candidate or edge_status == "candidate"
            if child not in reachable:
                reachable.add(child)
            child_state = (child, child_has_candidate)
            if child_state not in seen_states:
                seen_states.add(child_state)
                queue.append(child_state)
            if child_has_candidate:
                candidate_reachable.add(child)

    carbon_reactant_entities = set()
    for reaction in (e for e in network["entities"] if e.get("type") == "biochemical_reaction"):
        for entity_id in _participant_ids(network, reaction["id"], "has_reactant"):
            if _carbon_count(entities.get(entity_id, {}).get("attributes", {}).get("canonicalSmiles")):
                carbon_reactant_entities.add(entity_id)
    external_carbon_inputs = sorted(entity_id for entity_id in carbon_reactant_entities if entity_id != co2_id and not any(node[0] == entity_id for node in reachable))

    compounds_by_id = {c["id"]: c for c in compounds}
    exact = {row["terpedia_id"]: row for row in crosswalk["matches"]}
    targets = []
    for compound in compounds:
        match = next((row for row in crosswalk["matches"] if row["cannabisdb"]["cannabisdb_id"] == compound["id"]), None)
        if match is None:
            targets.append({"cannabisdb_id": compound["id"], "status": "unresolved", "reason": "no-exact-terpedia-identity", "carbon_atom_count": compound["carbon_atom_count"], "reachable_carbon_atoms": 0})
            continue
        entity_id = match["terpedia_id"]
        entity_smiles = entities.get(entity_id, {}).get("attributes", {}).get("canonicalSmiles")
        product_nodes = {(entity_id, i) for i in _carbon_atom_indices(entity_smiles)}
        reachable_count = len(product_nodes & reachable)
        if reachable_count == len(product_nodes) and reachable_count and not (product_nodes & candidate_reachable):
            status, reason = "supported", "all-entity-product-carbons-reachable-from-CO2"
        elif reachable_count:
            status, reason = "candidate", "partial-entity-product-carbon-reachability"
        else:
            status, reason = "unresolved", "entity-not-reachable-from-CO2-through-inferred-carbon-edges"
        targets.append({"cannabisdb_id": compound["id"], "terpedia_id": entity_id, "status": status, "reason": reason, "carbon_atom_count": compound["carbon_atom_count"], "entity_product_carbon_atoms": len(product_nodes), "reachable_carbon_atoms": reachable_count})
    report = {"schema": "cannabis-carbon.carbon-lineage.v1", "source": str(network_path), "direction_overrides": directions, "direction_conflicts": direction_conflicts, "carbon_source_policy": "CO2 is the only admissible carbon source for a Cannabis plant; every other carbon-containing reactant is an explicit external-carbon-source blocker until connected to CO2.", "co2_entity_id": co2_id, "resolved_carbon_edges": len(edges), "inferred_carbon_edges": sum(e["status"] == "inferred" for e in edges), "candidate_carbon_edges": sum(e["status"] == "candidate" for e in edges), "reachable_carbon_nodes": len(reachable), "reachable_carbon_entity_ids": sorted({node[0] for node in reachable}), "external_carbon_input_entity_count": len(external_carbon_inputs), "external_carbon_input_entity_ids": external_carbon_inputs, "edge_block_counts": dict(edge_blocks), "target_summary": {status: sum(t["status"] == status for t in targets) for status in ("supported", "candidate", "unresolved")}, "carbon_edges": edges, "targets": targets, "claim_boundary": "Reachability is based on one-to-one structural RDKit mappings and exact identity crosswalks. It is not isotope tracing, enzyme validation, or proof of in-vivo biosynthesis; all unresolved reasons are retained."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {"resolved_carbon_edges": len(edges), "inferred_carbon_edges": sum(e["status"] == "inferred" for e in edges), "candidate_carbon_edges": sum(e["status"] == "candidate" for e in edges), "reachable_carbon_nodes": len(reachable), "target_summary": report["target_summary"]}
