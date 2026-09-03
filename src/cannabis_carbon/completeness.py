from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


def compute_completeness(network_path: Path, compounds_path: Path, mapping_path: Path | None = None) -> dict:
    with __import__("gzip").open(network_path, "rt", encoding="utf-8") as handle:
        network = json.load(handle)
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
        "coverage": {"mapped_carbon_atoms": 0, "unresolved_carbon_atoms": sum(c["carbon_atom_count"] for c in compounds), "coverage_percent": 0.0},
        "claim_boundary": "These are database-coverage metrics, not evidence that every listed compound is biosynthesized by Cannabis.",
    }
    if mapping_path and mapping_path.exists():
        mapping = json.loads(mapping_path.read_text())
        mapped = sum(len([m for m in row["mappings"] if m["status"] == "inferred"]) for row in mapping["reactions"])
        unresolved = sum(len(row["unresolved_product_carbons"]) for row in mapping["reactions"])
        result["coverage"].update(mapped_carbon_atoms=mapped, unresolved_carbon_atoms=unresolved, reaction_mapping_status_counts=mapping["status_counts"])
    return result
