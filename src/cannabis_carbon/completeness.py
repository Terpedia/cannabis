from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from .terpedia import load_network


def compute_completeness(network_path: Path, compounds_path: Path, mapping_path: Path | None = None, crosswalk_path: Path | None = None, lineage_path: Path | None = None, atom_audit_path: Path | None = None) -> dict:
    network = load_network(network_path)
    entities = {e["id"]: e for e in network["entities"]}
    metabolites = {i for i, e in entities.items() if e.get("type") == "metabolite"}
    reactions = {i for i, e in entities.items() if e.get("type") == "biochemical_reaction"}
    reactant_metabolites, product_metabolites, enzyme_reactions = set(), set(), set()
    for s in network["statements"]:
        if s["predicate"] == "has_reactant" and s["subjectId"] in reactions: reactant_metabolites.add(s["objectEntityId"])
        elif s["predicate"] == "has_product" and s["subjectId"] in reactions: product_metabolites.add(s["objectEntityId"])
        elif s["predicate"] in ("catalyzes", "maps_to_reaction") and s["objectEntityId"] in reactions: enzyme_reactions.add(s["objectEntityId"])
    compounds = json.loads(compounds_path.read_text())["compounds"]
    result = {
        "schema": "cannabis-carbon.completeness.v1",
        "terpedia": {
            "metabolites_total": len(metabolites),
            "metabolites_in_reactions": len(reactant_metabolites | product_metabolites),
            "metabolites_without_reactions": len(metabolites - reactant_metabolites - product_metabolites),
            "reactions_total": len(reactions),
            "reactions_with_enzyme_association": len(enzyme_reactions),
            "reactions_without_enzyme_association": len(reactions - enzyme_reactions),
            "metabolite_ids_without_reactions": sorted(metabolites - reactant_metabolites - product_metabolites),
            "reaction_ids_without_enzyme_association": sorted(reactions - enzyme_reactions),
        },
        "cannabisdb": {
            "compounds_total": len(compounds),
            "carbon_atoms_total": sum(c["carbon_atom_count"] for c in compounds),
            "compound_to_terpedia_identity_crosswalk": "not_yet_available",
        },
        "coverage": {"mapped_carbon_atoms": 0, "unresolved_carbon_atoms": sum(c["carbon_atom_count"] for c in compounds), "coverage_percent": None, "coverage_denominator": "all CannabisDB carbons; no complete pathway crosswalk"},
        "claim_boundary": "These are database-coverage metrics, not evidence that every listed compound is biosynthesized by Cannabis.",
    }
    if crosswalk_path and crosswalk_path.exists():
        crosswalk = json.loads(crosswalk_path.read_text())
        matched_ids = {row["cannabisdb"]["cannabisdb_id"] for row in crosswalk["matches"]}
        matched_carbons = sum(c["carbon_atom_count"] for c in compounds if c["id"] in matched_ids)
        result["cannabisdb"].update(compounds_with_exact_terpedia_identity=len(matched_ids), compounds_with_ambiguous_identity=crosswalk["ambiguous"], compounds_without_exact_terpedia_identity=crosswalk.get("cannabisdb_unmatched", crosswalk["unmatched"]), connectivity_candidate_identity_links=crosswalk.get("connectivity_candidate_matches", 0), connectivity_candidate_ambiguous=crosswalk.get("connectivity_candidate_ambiguous", 0), crosswalk_matched_carbon_atoms=matched_carbons, compound_to_terpedia_identity_crosswalk="exact-inchikey")
    if mapping_path and mapping_path.exists():
        mapping = json.loads(mapping_path.read_text())
        mapped = mapping["carbon_counts"]["mapped_carbon_atoms"]
        unresolved = mapping["carbon_counts"]["unresolved_or_ambiguous_carbon_atoms"]
        result["coverage"].update(mapped_carbon_atoms=mapped, unresolved_carbon_atoms=unresolved, reaction_product_carbon_atoms=mapping["carbon_counts"]["product_carbon_atoms"], reaction_mapping_coverage_percent=mapping["carbon_counts"]["mapping_coverage_percent"], reaction_mapping_status_counts=mapping["status_counts"])
    if lineage_path and lineage_path.exists():
        lineage = json.loads(lineage_path.read_text())
        result["coverage"]["co2_lineage"] = {"target_summary": lineage["target_summary"], "reachable_carbon_nodes": lineage["reachable_carbon_nodes"], "resolved_carbon_edges": lineage["resolved_carbon_edges"], "inferred_carbon_edges": lineage["inferred_carbon_edges"], "candidate_carbon_edges": lineage["candidate_carbon_edges"], "external_carbon_input_entity_count": lineage["external_carbon_input_entity_count"], "carbon_source_policy": lineage["carbon_source_policy"]}
    if atom_audit_path and atom_audit_path.exists():
        atom_audit = json.loads(atom_audit_path.read_text())
        result["coverage"]["carbon_atom_audit"] = {"source": str(atom_audit_path), "carbon_atoms_total": atom_audit.get("carbon_atoms_total"), "status_counts": atom_audit.get("status_counts"), "compound_count": atom_audit.get("compound_count")}
    return result
