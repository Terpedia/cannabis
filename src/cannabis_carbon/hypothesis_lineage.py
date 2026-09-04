"""Candidate CO2 reachability through source-linked hypothesis edges."""

from __future__ import annotations

import json
from collections import Counter, deque
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def _carbon_count(smiles: str | None) -> int:
    molecule = Chem.MolFromSmiles(smiles or "")
    return sum(atom.GetAtomicNum() == 6 for atom in molecule.GetAtoms()) if molecule else 0


def _structure_key(smiles: str | None) -> str:
    molecule = Chem.MolFromSmiles(smiles or "")
    return Chem.MolToSmiles(molecule, canonical=True) if molecule else (smiles or "")


def _carbon_input_accounting(edge: dict, source_carbon_atoms: int, product_carbon_atoms: int) -> dict:
    """Expose carbon-bearing reaction inputs without claiming atom correspondence.

    The forward-connections table retains the full Rhea reaction SMARTS and the
    corpus substrates that were or were not resolved.  Counting those inputs
    lets the report distinguish a genuine carbon-source requirement from a
    mapping failure.  It deliberately does not infer which input atom becomes
    which product atom.
    """
    try:
        required = json.loads(edge.get("required_substrate_structures_json") or "[]")
    except (TypeError, json.JSONDecodeError):
        required = []
    try:
        missing = json.loads(edge.get("missing_corpus_substrates_json") or "[]")
    except (TypeError, json.JSONDecodeError):
        missing = []
    required_carbon_counts = [_carbon_count(smiles) for smiles in required]
    missing_carbon_counts = [_carbon_count(smiles) for smiles in missing]
    ancillary_required = sum(required_carbon_counts[1:]) if required_carbon_counts else 0
    missing_carbon = sum(missing_carbon_counts)
    missing_by_structure = Counter(_structure_key(smiles) for smiles in missing)
    input_sources = []
    for index, smiles in enumerate(required):
        carbon_count = required_carbon_counts[index]
        if not carbon_count:
            continue
        structure_key = _structure_key(smiles)
        missing_count = min(1, missing_by_structure.get(structure_key, 0))
        if missing_count:
            missing_by_structure[structure_key] -= missing_count
        input_sources.append({
            "input_index": index,
            "carbon_atom_count": carbon_count,
            "structure": smiles,
            "canonical_structure": structure_key,
            "corpus_status": "missing" if missing_count else "represented-or-unresolved-by-structure",
        })
    single_carbon_ancillary_inputs = [
        item for item in input_sources
        if item["input_index"] > 0 and item["carbon_atom_count"] == 1
    ]
    endpoint_delta = product_carbon_atoms - source_carbon_atoms
    if endpoint_delta > 0 and ancillary_required:
        status = "additional-carbon-input-required"
    elif endpoint_delta > 0:
        status = "product-carbon-deficit"
    elif missing_carbon:
        status = "missing-carbon-bearing-cosubstrate"
    else:
        status = "endpoint-carbon-count-compatible"
    return {
        "source_endpoint_carbon_atoms": source_carbon_atoms,
        "product_carbon_atoms": product_carbon_atoms,
        "endpoint_carbon_delta": endpoint_delta,
        "required_input_carbon_atoms": sum(required_carbon_counts),
        "required_ancillary_carbon_atoms": ancillary_required,
        "missing_input_carbon_atoms": missing_carbon,
        "required_input_carbon_counts": required_carbon_counts,
        "missing_input_carbon_counts": missing_carbon_counts,
        "carbon_bearing_input_sources": input_sources,
        "unique_single_carbon_ancillary_input": single_carbon_ancillary_inputs[0] if len(single_carbon_ancillary_inputs) == 1 else None,
        "status": status,
    }


def _carbon_mcs(source_smiles: str | None, product_smiles: str | None) -> dict:
    source = Chem.MolFromSmiles(source_smiles or "")
    product = Chem.MolFromSmiles(product_smiles or "")
    product_carbons = {a.GetIdx() for a in product.GetAtoms() if a.GetAtomicNum() == 6} if product else set()
    if source is None or product is None:
        return {"status": "unresolved", "product_carbon_atoms": len(product_carbons), "mapped_product_carbon_atoms": 0, "reason": "missing-endpoint-structure"}
    result = rdFMCS.FindMCS([source, product], atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareAny, ringMatchesRingOnly=False, completeRingsOnly=False, timeout=1)
    if result.canceled:
        return {"status": "unresolved", "product_carbon_atoms": len(product_carbons), "mapped_product_carbon_atoms": 0, "reason": "mcs-timeout"}
    query = Chem.MolFromSmarts(result.smartsString)
    product_matches = product.GetSubstructMatches(query, uniquify=True) if query else []
    mapped = {atom for match in product_matches for atom in match if atom in product_carbons}
    return {"status": "candidate" if mapped == product_carbons else "unresolved", "product_carbon_atoms": len(product_carbons), "mapped_product_carbon_atoms": len(mapped), "mcs_atoms": result.numAtoms, "mapping_alternatives": len(product_matches), "reason": "complete-product-carbon-coverage" if mapped == product_carbons else "partial-product-carbon-coverage"}


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
        atom_mapping = []
        for step in path:
            edge = next((candidate for candidate in candidate_edges if candidate.get("reaction_id") == step["reaction_id"] and candidate.get("substrate_compound_id") == step["from_compound_id"] and candidate.get("product_compound_id") == step["to_compound_id"]), None)
            if edge is None:
                atom_mapping.append({"reaction_id": step["reaction_id"], "status": "unresolved", "reason": "hypothesis-edge-not-found"})
                continue
            source_smiles = by_id.get(edge.get("substrate_compound_id"), {}).get("smiles")
            product_smiles = by_id.get(edge.get("product_compound_id"), {}).get("smiles")
            mapping = _carbon_mcs(source_smiles, product_smiles)
            mapping["carbon_input_accounting"] = _carbon_input_accounting(
                edge,
                mapping.get("source_carbon_atoms", _carbon_count(source_smiles)),
                mapping.get("product_carbon_atoms", _carbon_count(product_smiles)),
            )
            atom_mapping.append({"reaction_id": step["reaction_id"], **mapping})
        complete_atom_path = bool(path) and all(item.get("status") == "candidate" for item in atom_mapping)
        target_rows.append({
            "cannabisdb_id": compound_id,
            "label": compound.get("label"),
            "carbon_atom_count": compound.get("carbon_atom_count", 0),
            "core_reachable_carbon_atoms": core_atoms,
            "status": status,
            "reason": reason,
            "path": path,
            "atom_mapping": atom_mapping,
            "atom_mapping_status": "complete-candidate" if complete_atom_path else "not-complete",
        })

    status_counts = Counter(row["status"] for row in target_rows)
    carbon_counts = Counter()
    for row in target_rows:
        carbon_counts[row["status"]] += row["carbon_atom_count"]
    blocked_edges = Counter(edge.get("blocker") or "unspecified" for edge in edges if edge.get("status") == "unresolved")
    carbon_source_gaps = []
    seen_gaps = set()
    for target in target_rows:
        if target["status"] != "candidate":
            continue
        for step, mapping in zip(target["path"], target["atom_mapping"]):
            accounting = mapping.get("carbon_input_accounting")
            if not accounting:
                continue
            if accounting["status"] == "endpoint-carbon-count-compatible" and mapping.get("status") == "candidate":
                continue
            key = (target["cannabisdb_id"], step["reaction_id"], step["from_compound_id"], step["to_compound_id"])
            if key in seen_gaps:
                continue
            seen_gaps.add(key)
            carbon_source_gaps.append({
                "cannabisdb_id": target["cannabisdb_id"],
                "target_label": target["label"],
                "reaction_id": step["reaction_id"],
                "from_compound_id": step["from_compound_id"],
                "to_compound_id": step["to_compound_id"],
                "atom_mapping_status": mapping.get("status"),
                "atom_mapping_reason": mapping.get("reason"),
                "carbon_input_accounting": accounting,
                "source_url": step.get("source_url"),
            })
    gap_status_counts = Counter(row["carbon_input_accounting"]["status"] for row in carbon_source_gaps)
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
            "candidate_targets_with_complete_atom_paths": sum(row["atom_mapping_status"] == "complete-candidate" for row in target_rows),
            "candidate_carbon_atoms_with_complete_atom_paths": sum(row["carbon_atom_count"] for row in target_rows if row["atom_mapping_status"] == "complete-candidate"),
        },
        "blocked_unresolved_hypothesis_edges": dict(sorted(blocked_edges.items())),
        "carbon_source_gap_summary": {
            "path_steps_with_explicit_carbon_gaps": len(carbon_source_gaps),
            "status_counts": dict(sorted(gap_status_counts.items())),
            "missing_input_carbon_atoms": sum(row["carbon_input_accounting"]["missing_input_carbon_atoms"] for row in carbon_source_gaps),
            "product_carbon_deficit_atoms": sum(max(0, row["carbon_input_accounting"]["endpoint_carbon_delta"]) for row in carbon_source_gaps if row["carbon_input_accounting"]["status"] == "product-carbon-deficit"),
        },
        "carbon_source_gaps": carbon_source_gaps,
        "targets": target_rows,
        "claim_boundary": "Candidate paths are source-linked structural hypotheses, not confirmed enzyme activity, isotope tracing, or proof of in-vivo Cannabis biosynthesis. Unresolved edges are not traversed.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return report
