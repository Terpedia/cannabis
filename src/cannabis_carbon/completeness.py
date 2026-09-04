from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

from .terpedia import load_network


_SPECIALTY_NAME = re.compile(r"cannab|tetrahydrocannabin|cannabidiol|cannabiger|cannabichrom|cannabinol|cannabicycl|cannabielso|cannabifuran|cannabitriol|cannabid|cannabivarin|cannabistilbene|cannabisativine", re.IGNORECASE)


def compute_completeness(network_path: Path, compounds_path: Path, mapping_path: Path | None = None, crosswalk_path: Path | None = None, lineage_path: Path | None = None, atom_audit_path: Path | None = None, hypotheses_path: Path | None = None, pubchem_path: Path | None = None, networkdb_path: Path | None = None, hypothesis_lineage_path: Path | None = None, candidate_path_carbon_path: Path | None = None) -> dict:
    network = load_network(network_path)
    entities = {e["id"]: e for e in network["entities"]}
    metabolites = {i for i, e in entities.items() if e.get("type") == "metabolite"}
    reactions = {i for i, e in entities.items() if e.get("type") == "biochemical_reaction"}
    non_enzymatic_reactions = {i for i in reactions if str(entities[i].get("attributes", {}).get("reactionClass", "")).startswith("non-enzymatic-")}
    enzyme_requiring_reactions = reactions - non_enzymatic_reactions
    reactant_metabolites, product_metabolites, enzyme_reactions = set(), set(), set()
    for s in network["statements"]:
        if s["predicate"] == "has_reactant" and s["subjectId"] in reactions: reactant_metabolites.add(s["objectEntityId"])
        elif s["predicate"] == "has_product" and s["subjectId"] in reactions: product_metabolites.add(s["objectEntityId"])
        elif s["predicate"] in ("catalyzes", "maps_to_reaction", "has_catalytic_activity") and s["objectEntityId"] in reactions: enzyme_reactions.add(s["objectEntityId"])
    compounds = json.loads(compounds_path.read_text())["compounds"]
    candidate_reactions = set()
    if hypotheses_path and hypotheses_path.exists():
        hypotheses = json.loads(hypotheses_path.read_text())
        hypothesis_items = hypotheses.get("items") or hypotheses.get("hypotheses", [])
        candidate_reactions = {item.get("reaction_id") for item in hypothesis_items if item.get("reaction_id") and item.get("candidate_proteins")}
    no_enzyme_reactions = enzyme_requiring_reactions - enzyme_reactions
    result = {
        "schema": "cannabis-carbon.completeness.v1",
        "terpedia": {
            "metabolites_total": len(metabolites),
            "metabolites_in_reactions": len(reactant_metabolites | product_metabolites),
            "metabolites_without_reactions": len(metabolites - reactant_metabolites - product_metabolites),
            "reactions_total": len(reactions),
            "non_enzymatic_reactions": len(non_enzymatic_reactions),
            "enzyme_requiring_reactions": len(enzyme_requiring_reactions),
            "reactions_with_enzyme_association": len(enzyme_reactions),
            "reactions_without_enzyme_association": len(no_enzyme_reactions),
            "reactions_without_enzyme_with_candidate_proteins": len(no_enzyme_reactions & candidate_reactions),
            "reactions_without_enzyme_without_candidate_proteins": len(no_enzyme_reactions - candidate_reactions),
            "metabolite_ids_without_reactions": sorted(metabolites - reactant_metabolites - product_metabolites),
            "reaction_ids_without_enzyme_association": sorted(no_enzyme_reactions),
            "non_enzymatic_reaction_ids": sorted(non_enzymatic_reactions),
        },
        "cannabisdb": {
            "compounds_total": len(compounds),
            "carbon_atoms_total": sum(c["carbon_atom_count"] for c in compounds),
            "compound_to_terpedia_identity_crosswalk": "not_yet_available",
            "pubchem_resolution": "not_yet_available",
            "external_id_coverage": {"records_with_any_external_id": sum(bool(c.get("external_ids")) for c in compounds), "records_without_external_ids": sum(not c.get("external_ids") for c in compounds), "counts_by_database": dict(sorted({db: sum(db in c.get("external_ids", {}) for c in compounds) for db in {db for c in compounds for db in c.get("external_ids", {})}}.items()))},
        },
        "coverage": {"mapped_carbon_atoms": 0, "unresolved_carbon_atoms": sum(c["carbon_atom_count"] for c in compounds), "coverage_percent": None, "coverage_denominator": "all CannabisDB carbons; no complete pathway crosswalk"},
        "claim_boundary": "These are database-coverage metrics, not evidence that every listed compound is biosynthesized by Cannabis.",
    }
    if crosswalk_path and crosswalk_path.exists():
        crosswalk = json.loads(crosswalk_path.read_text())
        matched_ids = {row["cannabisdb"]["cannabisdb_id"] for row in crosswalk["matches"]}
        matched_carbons = sum(c["carbon_atom_count"] for c in compounds if c["id"] in matched_ids)
        result["cannabisdb"].update(compounds_with_exact_terpedia_identity=len(matched_ids), compounds_with_ambiguous_identity=crosswalk["ambiguous"], compounds_without_exact_terpedia_identity=len(compounds) - len(matched_ids), compounds_without_any_identity_resolution=crosswalk.get("cannabisdb_unmatched", crosswalk["unmatched"]), connectivity_candidate_identity_links=crosswalk.get("connectivity_candidate_matches", 0), connectivity_candidate_ambiguous=crosswalk.get("connectivity_candidate_ambiguous", 0), tautomer_candidate_identity_links=crosswalk.get("tautomer_candidate_matches", 0), tautomer_candidate_ambiguous=crosswalk.get("tautomer_candidate_ambiguous", 0), name_candidate_identity_links=crosswalk.get("name_candidate_matches", 0), name_candidate_ambiguous=crosswalk.get("name_candidate_ambiguous", 0), crosswalk_matched_carbon_atoms=matched_carbons, compound_to_terpedia_identity_crosswalk="exact-inchikey-plus-candidate-tautomer-and-formula-compatible-name-layers")
    if pubchem_path and pubchem_path.exists():
        pubchem = json.loads(pubchem_path.read_text())
        result["cannabisdb"]["pubchem_resolution"] = pubchem.get("summary", {})
    if mapping_path and mapping_path.exists():
        mapping = json.loads(mapping_path.read_text())
        mapped = mapping["carbon_counts"]["mapped_carbon_atoms"]
        unresolved = mapping["carbon_counts"]["unresolved_or_ambiguous_carbon_atoms"]
        mapping_rows = [mapping_row for reaction in mapping.get("reactions", []) for mapping_row in reaction.get("mappings", [])]
        row_status_counts = Counter(row.get("status", "unresolved") for row in mapping_rows)
        blocked_reactions = [reaction for reaction in mapping.get("reactions", []) if any(row.get("status") in ("ambiguous", "unresolved") for row in reaction.get("mappings", []))]
        result["coverage"].update(mapped_carbon_atoms=mapped, unresolved_carbon_atoms=unresolved, reaction_product_carbon_atoms=mapping["carbon_counts"]["product_carbon_atoms"], reaction_mapping_coverage_percent=mapping["carbon_counts"]["mapping_coverage_percent"], reaction_mapping_status_counts=mapping["status_counts"], carbon_mapping_blockers={"reactions_with_blocked_product_carbon_rows": len(blocked_reactions), "product_carbon_row_status_counts": dict(sorted(row_status_counts.items())), "blocked_product_carbon_rows": sum(row_status_counts.get(status, 0) for status in ("ambiguous", "unresolved"))})
    if lineage_path and lineage_path.exists():
        lineage = json.loads(lineage_path.read_text())
        status_identity = Counter()
        carbon_status = Counter()
        specialty_status = Counter()
        blocking_reasons = Counter()
        compound_index = {c.get("id"): c for c in compounds}
        for target in lineage.get("targets", []):
            status = target.get("status", "unresolved")
            identity = target.get("identity_status") or "unresolved"
            status_identity[f"{status}:{identity}"] += 1
            carbon_status[status] += target.get("carbon_atom_count", 0)
            record = compound_index.get(target.get("cannabisdb_id"), {})
            if any(_SPECIALTY_NAME.search(name or "") for name in [record.get("label", ""), *record.get("aliases", [])]):
                specialty_status[status] += 1
            blocking_reasons[target.get("reason") or "unspecified"] += 1
        result["coverage"]["co2_lineage"] = {"target_summary": lineage["target_summary"], "reachable_carbon_nodes": lineage["reachable_carbon_nodes"], "resolved_carbon_edges": lineage["resolved_carbon_edges"], "inferred_carbon_edges": lineage["inferred_carbon_edges"], "candidate_carbon_edges": lineage["candidate_carbon_edges"], "external_carbon_input_entity_count": lineage["external_carbon_input_entity_count"], "carbon_source_policy": lineage["carbon_source_policy"], "target_triage": {"target_counts_by_status_and_identity": dict(sorted(status_identity.items())), "carbon_atoms_by_target_status": dict(sorted(carbon_status.items())), "specialty_target_counts_by_status": dict(sorted(specialty_status.items())), "blocking_reason_counts": dict(sorted(blocking_reasons.items()))}}
    if atom_audit_path and atom_audit_path.exists():
        atom_audit = json.loads(atom_audit_path.read_text())
        audit_counts = atom_audit.get("status_counts", {})
        audit_total = atom_audit.get("carbon_atoms_total") or 0
        evidence_atoms = sum(audit_counts.get(status, 0) for status in ("supported", "candidate", "inferred"))
        result["coverage"]["carbon_atom_audit"] = {"source": str(atom_audit_path), "carbon_atoms_total": audit_total, "status_counts": audit_counts, "compound_count": atom_audit.get("compound_count"), "strict_supported_percent": round(100 * audit_counts.get("supported", 0) / audit_total, 6) if audit_total else None, "evidence_bearing_percent": round(100 * evidence_atoms / audit_total, 6) if audit_total else None, "unresolved_percent": round(100 * audit_counts.get("unresolved", 0) / audit_total, 6) if audit_total else None}
    if networkdb_path and networkdb_path.exists():
        networkdb = json.loads(networkdb_path.read_text())
        hypothesis_edges = networkdb.get("hypothetical_connections", [])
        result["hypothesis_layer"] = {
            "source": str(networkdb_path),
            "compound_records": len(networkdb.get("compounds", [])),
            "hypothesis_edges": len(hypothesis_edges),
            "hypothesis_edge_status_counts": dict(sorted(Counter(edge.get("status", "unresolved") for edge in hypothesis_edges).items())),
            "hypothesis_edges_with_enzyme_evidence": sum(bool(edge.get("enzyme_evidence") or edge.get("enzyme_catalog")) for edge in hypothesis_edges),
            "hypothesis_reaction_inventory": networkdb.get("coverage", {}).get("hypothetical_reaction_inventory", 0),
            "hypothesis_product_inventory": networkdb.get("coverage", {}).get("hypothetical_product_inventory", 0),
            "missing_substrate_nodes": networkdb.get("coverage", {}).get("hypothetical_missing_substrate_nodes", 0),
            "claim_boundary": "Hypothesis-layer counts are separate from balanced reaction and CO2 lineage metrics; they represent source-linked candidates and unresolved structures, not confirmed Cannabis metabolism.",
        }
    if hypothesis_lineage_path and hypothesis_lineage_path.exists():
        hypothesis_lineage = json.loads(hypothesis_lineage_path.read_text())
        result.setdefault("hypothesis_layer", {})["candidate_lineage"] = {
            "source": str(hypothesis_lineage_path),
            "seed_core_reachable_compounds": hypothesis_lineage.get("seed_core_reachable_compounds", 0),
            "candidate_edges_traversed": hypothesis_lineage.get("candidate_edges_traversed", 0),
            "candidate_reachable_entities": hypothesis_lineage.get("candidate_reachable_entities", 0),
            "target_summary": hypothesis_lineage.get("target_summary", {}),
            "blocked_unresolved_hypothesis_edges": hypothesis_lineage.get("blocked_unresolved_hypothesis_edges", {}),
            "claim_boundary": hypothesis_lineage.get("claim_boundary"),
        }
    if candidate_path_carbon_path and candidate_path_carbon_path.exists():
        candidate_paths = json.loads(candidate_path_carbon_path.read_text())
        path_rows = candidate_paths.get("rows", [])
        candidate_statuses = Counter(row.get("carbon_mapping", {}).get("status", "unresolved") for row in path_rows)
        candidate_products = {row.get("candidate_product_terpene_id") for row in path_rows if row.get("candidate_product_terpene_id")}
        candidate_core_edges = sum(len(row.get("core_path_carbon_edges", [])) for row in path_rows)
        result["coverage"]["candidate_co2_path_carbon_layer"] = {
            "source": str(candidate_path_carbon_path),
            "path_count": candidate_paths.get("path_count", len(path_rows)),
            "candidate_product_count": len(candidate_products),
            "complete_product_carbon_paths": candidate_paths.get("paths_with_complete_product_carbon_mapping", 0),
            "mapped_product_carbon_atoms": candidate_paths.get("mapped_product_carbon_atoms", 0),
            "unresolved_product_carbon_atoms": candidate_paths.get("unresolved_product_carbon_atoms", 0),
            "carbon_mapping_status_counts": dict(sorted(candidate_statuses.items())),
            "atom_level_core_path_edges": candidate_core_edges,
            "path_mode": "all-reactions-reversible-upper-bound",
            "claim_boundary": "These are source-linked, structurally mapped candidate paths. Reversible traversal is an upper bound and does not establish physiological direction, enzyme activity, isotope tracing, or in-vivo Cannabis biosynthesis.",
        }
    return result
