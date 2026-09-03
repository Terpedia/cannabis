"""Build a unified, source-preserving NetworkDB snapshot."""

from __future__ import annotations

import json
from pathlib import Path

from rdkit import Chem

from .terpedia import load_network


def build_networkdb(network_path: Path, compounds_path: Path, crosswalk_path: Path, output: Path, hypotheses_path: Path | None = None) -> dict:
    network = load_network(network_path)
    catalog = json.loads(compounds_path.read_text())["compounds"]
    crosswalk = json.loads(crosswalk_path.read_text())
    hypotheses = json.loads(hypotheses_path.read_text()) if hypotheses_path and hypotheses_path.exists() else {"items": []}
    entities = {e["id"]: e for e in network["entities"]}
    identity_by_cdb = {row["cannabisdb"]["cannabisdb_id"]: row["terpedia_id"] for row in crosswalk["matches"]}
    compounds = []
    for compound in catalog:
        compounds.append({**compound, "namespace": "cannabisdb", "identity_link": identity_by_cdb.get(compound["id"])})
    for entity in network["entities"]:
        if entity.get("type") != "metabolite":
            continue
        attrs = entity.get("attributes", {})
        mol = Chem.MolFromSmiles(attrs["canonicalSmiles"]) if attrs.get("canonicalSmiles") else None
        compounds.append({"id": entity["id"], "namespace": "terpedia", "label": entity.get("label", entity["id"]), "formula": attrs.get("molecularFormula"), "smiles": attrs.get("canonicalSmiles"), "source_url": entity.get("url"), "carbon_atom_count": sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()) if mol is not None else None})
    candidate_by_reaction = {}
    for item in hypotheses.get("items", []):
        if item.get("reaction_id"):
            candidate_by_reaction.setdefault(item["reaction_id"], []).extend(item.get("candidate_proteins", []))
    reactions = []
    for entity in network["entities"]:
        if entity.get("type") != "biochemical_reaction":
            continue
        reaction_id = entity["id"]
        statements = [s for s in network["statements"] if s.get("subjectId") == reaction_id]
        participants = lambda predicate: [{"compound_id": s["objectEntityId"], "coefficient": (s.get("qualifiers") or {}).get("stoichiometricCoefficient", 1), "compartment": (s.get("qualifiers") or {}).get("compartment")} for s in statements if s.get("predicate") == predicate]
        enzyme_statements = [s for s in statements if s.get("predicate") in ("catalyzes", "maps_to_reaction")]
        enzymes = sorted({s.get("subjectId") for s in enzyme_statements})
        attrs = entity.get("attributes", {})
        candidate_proteins = candidate_by_reaction.get(reaction_id, [])
        reactions.append({"id": reaction_id, "label": entity.get("label", reaction_id), "equation": attrs.get("equation"), "reaction_smiles": attrs.get("reactionSmiles"), "ec_numbers": attrs.get("ecNumbers", []), "reactants": participants("has_reactant"), "products": participants("has_product"), "enzyme_ids": enzymes, "enzyme_associations": [{"enzyme_id": s.get("subjectId"), "predicate": s.get("predicate"), "sources": s.get("sources", []), "qualifiers": s.get("qualifiers", {})} for s in enzyme_statements], "candidate_proteins": candidate_proteins, "status": "supported" if any((s.get("qualifiers") or {}).get("directExperimentalEvidence") for s in enzyme_statements) else "candidate" if enzymes or candidate_proteins else "unresolved", "source_url": entity.get("url"), "directional_rhea_ids": entity.get("identifiers", {}).get("directionalRheaIds", [])})
    candidate_protein_records = {p.get("proteinId"): p for item in hypotheses.get("items", []) for p in item.get("candidate_proteins", []) if p.get("proteinId")}
    report = {"schema": "cannabis-carbon.networkdb.v1", "sources": {"cannabisdb_compounds": str(compounds_path), "terpedia_network": str(network_path), "identity_crosswalk": str(crosswalk_path), "candidate_hypotheses": str(hypotheses_path) if hypotheses_path else None}, "carbon_source_policy": "CO2 is the only admissible carbon source for Cannabis; no carbon-containing compound is treated as an implicit source.", "compounds": compounds, "reactions": reactions, "identity_links": crosswalk["matches"], "candidate_proteins": list(candidate_protein_records.values()), "candidate_hypotheses": hypotheses.get("items", []), "coverage": {"cannabisdb_compounds": len(catalog), "terpedia_metabolites": sum(e.get("type") == "metabolite" for e in network["entities"]), "compound_records": len(compounds), "terpedia_reactions": len(reactions), "reaction_records": len(reactions), "exact_identity_links": len(crosswalk["matches"]), "ambiguous_identity_links": crosswalk["ambiguous"], "unmatched_cannabisdb_compounds": crosswalk["unmatched"], "candidate_hypotheses": len(hypotheses.get("items", [])), "candidate_protein_records": len(candidate_protein_records), "reactions_with_candidate_proteins": sum(bool(r["candidate_proteins"]) for r in reactions)}, "claim_boundary": "This is a unified inventory and reaction database. Presence of a compound, reaction, enzyme association, or identity link does not establish in-vivo cannabis biosynthesis; candidate proteins remain hypotheses."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return report["coverage"]
