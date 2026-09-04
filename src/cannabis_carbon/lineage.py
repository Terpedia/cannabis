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


def _entity_atom_index_map(reaction_molecule, entity_smiles: str | None) -> dict[int, int] | None:
    """Map reaction-SMILES atom indices onto the resolved entity structure.

    Rhea reaction SMILES and Terpedia canonical SMILES describe the same
    molecule with independent RDKit atom orderings.  A carbon edge must use a
    stable entity-local atom index, otherwise edges from different reactions
    can be joined to the wrong carbon.  Full-structure matching is required;
    failure remains an explicit participant-structure blocker.
    """
    entity_molecule = Chem.MolFromSmiles(entity_smiles) if entity_smiles else None
    if entity_molecule is None:
        return None
    match = entity_molecule.GetSubstructMatch(reaction_molecule)
    if len(match) != reaction_molecule.GetNumAtoms():
        return None
    return {reaction_index: entity_index for reaction_index, entity_index in enumerate(match)}


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
        reactant_atom_maps = [_entity_atom_index_map(molecule, entities.get(entity_id, {}).get("attributes", {}).get("canonicalSmiles")) if entity_id else None for molecule, entity_id in zip(reactants, reactant_ids)]
        product_atom_maps = [_entity_atom_index_map(molecule, entities.get(entity_id, {}).get("attributes", {}).get("canonicalSmiles")) if entity_id else None for molecule, entity_id in zip(products, product_ids)]
        mappings = map_reaction_smiles(reaction_smiles)["mappings"] if direction.get("orientation") == "reverse_master" else row["mappings"]
        for mapping in mappings:
            if mapping.get("status") not in ("inferred", "candidate", "ambiguous"):
                continue
            alternatives = mapping.get("alternatives") or [{"reactant_index": mapping.get("reactant_index"), "reactant_atom": mapping.get("reactant_atom")}]
            for alternative in alternatives:
                ri, pi = alternative["reactant_index"], mapping["product_index"]
                if ri is None or ri >= len(reactant_ids) or pi >= len(product_ids) or not reactant_ids[ri] or not product_ids[pi] or not reactant_atom_maps[ri] or not product_atom_maps[pi] or alternative["reactant_atom"] not in reactant_atom_maps[ri] or mapping["product_atom"] not in product_atom_maps[pi]:
                    edge_blocks["participant-structure-unresolved"] += 1
                    continue
                edges.append({"reaction_id": reaction_id, "reactant_entity_id": reactant_ids[ri], "reactant_atom": reactant_atom_maps[ri][alternative["reactant_atom"]], "product_entity_id": product_ids[pi], "product_atom": product_atom_maps[pi][mapping["product_atom"]], "status": "candidate" if mapping["status"] == "ambiguous" else mapping["status"], "provenance": direction.get("source") or row.get("rhea_url"), "directional_rhea_id": direction.get("directional_rhea_id"), "mapping_reason": mapping.get("reason")})

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
    # CO2 is the explicit seed.  Decarboxylation reactions can also emit CO2
    # through candidate edges, but that downstream occurrence must not demote
    # the seed's own identity from supported to candidate.
    candidate_reachable.difference_update(co2_atoms)

    carbon_reactant_entities = set()
    for reaction in (e for e in network["entities"] if e.get("type") == "biochemical_reaction"):
        for entity_id in _participant_ids(network, reaction["id"], "has_reactant"):
            if _carbon_count(entities.get(entity_id, {}).get("attributes", {}).get("canonicalSmiles")):
                carbon_reactant_entities.add(entity_id)
    external_carbon_inputs = sorted(entity_id for entity_id in carbon_reactant_entities if entity_id != co2_id and not any(node[0] == entity_id for node in reachable))

    exact_by_cdb = {row["cannabisdb"]["cannabisdb_id"]: row for row in crosswalk["matches"]}
    candidate_by_cdb = defaultdict(list)
    for row in crosswalk.get("candidate_matches", []):
        candidate_by_cdb[row["cannabisdb"]["cannabisdb_id"]].append(row)
    targets = []
    for compound in compounds:
        match = exact_by_cdb.get(compound["id"])
        if match is None:
            identity_candidates = candidate_by_cdb.get(compound["id"], [])
            if len(identity_candidates) != 1:
                targets.append({"cannabisdb_id": compound["id"], "status": "unresolved", "reason": "ambiguous-connectivity-identity" if len(identity_candidates) > 1 else "no-terpedia-identity", "identity_candidates": [{"terpedia_id": row["terpedia_id"], "label": row.get("terpedia_label"), "method": row.get("method")} for row in identity_candidates], "carbon_atom_count": compound["carbon_atom_count"], "reachable_carbon_atoms": 0})
                continue
            match = identity_candidates[0]
            identity_status = "candidate"
        else:
            identity_status = "exact"
        entity_id = match["terpedia_id"]
        entity_smiles = entities.get(entity_id, {}).get("attributes", {}).get("canonicalSmiles")
        product_nodes = {(entity_id, i) for i in _carbon_atom_indices(entity_smiles)}
        reachable_count = len(product_nodes & reachable)
        if identity_status == "candidate":
            status = "candidate" if reachable_count else "unresolved"
            reason = "candidate-connectivity-identity-and-all-entity-product-carbons-reachable-from-CO2" if reachable_count == len(product_nodes) and reachable_count else "candidate-connectivity-identity-with-partial-co2-lineage" if reachable_count else "candidate-connectivity-identity-not-reachable-from-co2"
        elif reachable_count == len(product_nodes) and reachable_count and not (product_nodes & candidate_reachable):
            status, reason = "supported", "all-entity-product-carbons-reachable-from-CO2"
        elif reachable_count:
            status, reason = "candidate", "partial-entity-product-carbon-reachability"
        else:
            status, reason = "unresolved", "entity-not-reachable-from-CO2-through-inferred-carbon-edges"
        targets.append({"cannabisdb_id": compound["id"], "terpedia_id": entity_id, "identity_status": identity_status, "status": status, "reason": reason, "carbon_atom_count": compound["carbon_atom_count"], "entity_product_carbon_atoms": len(product_nodes), "reachable_carbon_atoms": reachable_count})
    report = {"schema": "cannabis-carbon.carbon-lineage.v1", "source": str(network_path), "direction_overrides": directions, "direction_conflicts": direction_conflicts, "carbon_source_policy": "CO2 is the only admissible carbon source for a Cannabis plant; every other carbon-containing reactant is an explicit external-carbon-source blocker until connected to CO2.", "co2_entity_id": co2_id, "resolved_carbon_edges": len(edges), "inferred_carbon_edges": sum(e["status"] == "inferred" for e in edges), "candidate_carbon_edges": sum(e["status"] == "candidate" for e in edges), "reachable_carbon_nodes": len(reachable), "reachable_carbon_entity_ids": sorted({node[0] for node in reachable}), "external_carbon_input_entity_count": len(external_carbon_inputs), "external_carbon_input_entity_ids": external_carbon_inputs, "edge_block_counts": dict(edge_blocks), "target_summary": {status: sum(t["status"] == status for t in targets) for status in ("supported", "candidate", "unresolved")}, "carbon_edges": edges, "targets": targets, "claim_boundary": "Reachability is based on one-to-one structural RDKit mappings and exact identity crosswalks. It is not isotope tracing, enzyme validation, or proof of in-vivo biosynthesis; all unresolved reasons are retained."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {"resolved_carbon_edges": len(edges), "inferred_carbon_edges": sum(e["status"] == "inferred" for e in edges), "candidate_carbon_edges": sum(e["status"] == "candidate" for e in edges), "reachable_carbon_nodes": len(reachable), "target_summary": report["target_summary"]}


def build_carbon_atom_audit(network_path: Path, lineage_path: Path, crosswalk_path: Path, compounds_path: Path, output: Path) -> dict:
    """Partition every CannabisDB carbon atom into an auditable status group.

    The lineage report stores reaction-level edges; this companion artifact
    projects those edges back onto each CannabisDB structure after explicitly
    resolving atom-order differences.  Groups contain atom-index lists whose
    union is exactly the molecule's carbon atom set, so omitted atoms are
    detectable by validation rather than silently disappearing.
    """
    network = load_network(network_path)
    lineage = json.loads(lineage_path.read_text())
    crosswalk = json.loads(crosswalk_path.read_text())
    compounds = json.loads(compounds_path.read_text())["compounds"]
    entities = {e["id"]: e for e in network["entities"]}
    co2_id = lineage.get("co2_entity_id", "chebi:16526")
    co2_smiles = entities.get(co2_id, {}).get("attributes", {}).get("canonicalSmiles")
    co2_nodes = {(co2_id, i) for i in _carbon_atom_indices(co2_smiles)}
    reachable = set(co2_nodes)
    forward = defaultdict(list)
    for edge in lineage.get("carbon_edges", []):
        forward[(edge["reactant_entity_id"], edge["reactant_atom"])].append(edge)
    candidate_reachable = set()
    queue = deque((node, False) for node in co2_nodes)
    seen_states = {(node, False) for node in co2_nodes}
    while queue:
        node, has_candidate = queue.popleft()
        for edge in forward.get(node, []):
            child = (edge["product_entity_id"], edge["product_atom"])
            child_has_candidate = has_candidate or edge.get("status") == "candidate"
            reachable.add(child)
            if child_has_candidate:
                candidate_reachable.add(child)
            state = (child, child_has_candidate)
            if state not in seen_states:
                seen_states.add(state)
                queue.append(state)
    candidate_reachable.difference_update(co2_nodes)
    incoming = defaultdict(list)
    for edge in lineage.get("carbon_edges", []):
        node = (edge["product_entity_id"], edge["product_atom"])
        incoming[node].append(edge)

    exact_by_cdb = {row["cannabisdb"]["cannabisdb_id"]: row for row in crosswalk.get("matches", [])}
    candidate_by_cdb = defaultdict(list)
    for row in crosswalk.get("candidate_matches", []):
        candidate_by_cdb[row["cannabisdb"]["cannabisdb_id"]].append(row)

    def source_urls(*values):
        return sorted({str(value) for value in values if value})

    atom_groups = []
    status_counts = defaultdict(int)
    for compound in compounds:
        compound_molecule = Chem.MolFromSmiles(compound.get("smiles", ""))
        carbon_indices = [atom.GetIdx() for atom in compound_molecule.GetAtoms() if atom.GetAtomicNum() == 6] if compound_molecule else []
        match = exact_by_cdb.get(compound["id"])
        identity_status = "exact"
        if match is None:
            candidates = candidate_by_cdb.get(compound["id"], [])
            if len(candidates) == 1:
                match = candidates[0]
                identity_status = "candidate"
            else:
                match = None
                identity_status = None
        terpedia_id = match.get("terpedia_id") if match else None
        entity_smiles = entities.get(terpedia_id, {}).get("attributes", {}).get("canonicalSmiles") if terpedia_id else None
        index_map = _entity_atom_index_map(compound_molecule, entity_smiles) if compound_molecule is not None and entity_smiles else None
        grouped = {}
        for atom_index in carbon_indices:
            status = "unresolved"
            reason = "no-terpedia-identity"
            entity_atom = None
            reaction_ids = []
            provenance = source_urls(compound.get("source_url"), crosswalk_path)
            if terpedia_id and index_map is None:
                reason = "identity-atom-order-unresolved"
                provenance = source_urls(compound.get("source_url"), match.get("method"), entities.get(terpedia_id, {}).get("url"))
            elif terpedia_id:
                entity_atom = index_map.get(atom_index)
                node = (terpedia_id, entity_atom) if entity_atom is not None else None
                edges = incoming.get(node, []) if node else []
                reaction_ids = sorted({edge["reaction_id"] for edge in edges})
                provenance = source_urls(*(edge.get("provenance") for edge in edges), entities.get(terpedia_id, {}).get("url"), compound.get("source_url"))
                if node in co2_nodes:
                    status, reason = "supported", "co2-seed"
                elif node not in reachable:
                    reason = "entity-carbon-not-reachable-from-co2"
                elif identity_status == "candidate" or node in candidate_reachable:
                    status, reason = "candidate", "candidate-co2-lineage-or-identity"
                else:
                    status, reason = "inferred", "rdkit-structural-co2-lineage"
            key = (status, reason, tuple(reaction_ids), tuple(provenance), entity_atom)
            grouped.setdefault(key, {"status": status, "reason": reason, "atom_indices": [], "entity_atom_indices": [], "reaction_ids": reaction_ids, "provenance": provenance})
            grouped[key]["atom_indices"].append(atom_index)
            grouped[key]["entity_atom_indices"].append(entity_atom)
            status_counts[status] += 1
        atom_groups.append({"cannabisdb_id": compound["id"], "atom_index_namespace": "RDKit atom indices in the CannabisDB SMILES field", "carbon_atom_count": len(carbon_indices), "identity_status": identity_status, "terpedia_id": terpedia_id, "groups": list(grouped.values()), "claim_boundary": "Each group explicitly accounts for listed CannabisDB carbon atom indices; unresolved atoms are not assigned an origin."})
    report = {"schema": "cannabis-carbon.carbon-atom-audit.v1", "source_network": str(network_path), "source_lineage": str(lineage_path), "source_crosswalk": str(crosswalk_path), "source_compounds": str(compounds_path), "atom_index_namespace": "RDKit atom indices in each compound's CannabisDB SMILES field; these are not assumed to equal source-SDF atom ordering", "carbon_source_policy": lineage.get("carbon_source_policy"), "compound_count": len(compounds), "carbon_atoms_total": sum(item["carbon_atom_count"] for item in atom_groups), "status_counts": {status: status_counts[status] for status in ("supported", "candidate", "inferred", "unresolved")}, "compounds": atom_groups, "claim_boundary": "This is a structure-indexed provenance audit. Inferred and candidate statuses are not isotope tracing, enzyme validation, or proof of in-vivo cannabis biosynthesis."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {"compound_count": report["compound_count"], "carbon_atoms_total": report["carbon_atoms_total"], "status_counts": report["status_counts"]}
