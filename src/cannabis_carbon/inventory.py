"""Build named CannabisDB specialty-metabolite inventories without identity guessing."""

from __future__ import annotations

import json
import re
from pathlib import Path

from .terpedia import load_network


SPECIALTY_NAME = re.compile(
    r"cannab|tetrahydrocannabin|cannabidiol|cannabiger|cannabichrom|"
    r"cannabinol|cannabicycl|cannabielso|cannabifuran|cannabitriol|"
    r"cannabid|cannabivarin|cannabit|cannabisp|cannabistilbene|cannabisativine",
    re.IGNORECASE,
)


def build_specialty_inventory(compounds_path: Path, crosswalk_path: Path, network_path: Path, output: Path) -> dict:
    compounds = json.loads(compounds_path.read_text())["compounds"]
    crosswalk = json.loads(crosswalk_path.read_text())
    network = load_network(network_path)
    exact_by_cdb = {row["cannabisdb"]["cannabisdb_id"]: row for row in crosswalk.get("matches", [])}
    candidate_by_cdb = {}
    for row in crosswalk.get("candidate_matches", []):
        candidate_by_cdb.setdefault(row["cannabisdb"]["cannabisdb_id"], []).append(row)
    reactions_by_metabolite = {}
    for statement in network["statements"]:
        if statement.get("predicate") not in ("has_reactant", "has_product"):
            continue
        reactions_by_metabolite.setdefault(statement["objectEntityId"], []).append({
            "reaction_id": statement["subjectId"],
            "role": "reactant" if statement["predicate"] == "has_reactant" else "product",
            "coefficient": (statement.get("qualifiers") or {}).get("stoichiometricCoefficient", 1),
        })
    records = []
    for compound in compounds:
        names = [compound.get("label", ""), *compound.get("aliases", [])]
        if not any(SPECIALTY_NAME.search(name) for name in names):
            continue
        exact = exact_by_cdb.get(compound["id"])
        terpedia_id = exact["terpedia_id"] if exact else None
        reaction_entity_id = terpedia_id or f"cannabisdb:{compound['id']}"
        records.append({
            "cannabisdb_id": compound["id"],
            "label": compound.get("label", compound["id"]),
            "aliases": compound.get("aliases", []),
            "formula": compound.get("formula"),
            "smiles": compound.get("smiles"),
            "carbon_atom_count": compound.get("carbon_atom_count"),
            "identity_status": "exact" if exact else "candidate" if candidate_by_cdb.get(compound["id"]) else "unresolved",
            "terpedia_id": terpedia_id,
            "identity_candidates": [{"terpedia_id": row["terpedia_id"], "terpedia_label": row.get("terpedia_label")} for row in candidate_by_cdb.get(compound["id"], [])],
            "reaction_participation": reactions_by_metabolite.get(reaction_entity_id, []),
            "carbon_status": compound.get("carbon_status", {}),
            "source": compound.get("source"),
            "source_url": compound.get("source_url"),
        })
    report = {
        "schema": "cannabis-carbon.named-specialty-inventory.v1",
        "selection": "Source CannabisDB names matching the explicit cannabinoid/cannabis specialty-name regular expression; this is a review inventory, not a chemical-class assertion.",
        "source": str(compounds_path),
        "crosswalk_source": str(crosswalk_path),
        "network_source": str(network_path),
        "record_count": len(records),
        "carbon_atom_count": sum(record.get("carbon_atom_count") or 0 for record in records),
        "identity_status_counts": {status: sum(record["identity_status"] == status for record in records) for status in ("exact", "candidate", "unresolved")},
        "records_with_reaction_participation": sum(bool(record["reaction_participation"]) for record in records),
        "records": records,
        "claim_boundary": "Name matching identifies records for curation. It does not establish that a compound is endogenous, that a reaction produces it in Cannabis, or that any candidate protein is functional.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {key: report[key] for key in ("record_count", "carbon_atom_count", "identity_status_counts", "records_with_reaction_participation")}
