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


def build_specialty_inventory(
    compounds_path: Path,
    crosswalk_path: Path,
    network_path: Path,
    output: Path,
    lineage_path: Path | None = None,
    pubchem_path: Path | None = None,
) -> dict:
    compounds = json.loads(compounds_path.read_text())["compounds"]
    crosswalk = json.loads(crosswalk_path.read_text())
    network = load_network(network_path)
    lineage = json.loads(lineage_path.read_text()) if lineage_path and lineage_path.exists() else {}
    lineage_by_cdb = {row["cannabisdb_id"]: row for row in lineage.get("targets", [])}
    pubchem = json.loads(pubchem_path.read_text()) if pubchem_path and pubchem_path.exists() else {}
    pubchem_by_cdb = {row["cannabisdb_id"]: row for row in pubchem.get("records", [])}
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
        lineage_target = lineage_by_cdb.get(compound["id"], {})
        pubchem_record = pubchem_by_cdb.get(compound["id"], {})
        terpedia_id = exact["terpedia_id"] if exact else None
        reaction_entity_id = terpedia_id or f"cannabisdb:{compound['id']}"
        unresolved_carbon_atoms = max(
            (compound.get("carbon_atom_count") or 0) - (lineage_target.get("reachable_carbon_atoms") or 0),
            0,
        )
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
            "carbon_lineage_status": lineage_target.get("status", "unresolved"),
            "reachable_carbon_atoms": lineage_target.get("reachable_carbon_atoms", 0),
            "unresolved_carbon_atoms": unresolved_carbon_atoms,
            "reversible_upper_bound_reachable_carbon_atoms": lineage_target.get(
                "reversible_upper_bound_reachable_carbon_atoms", 0
            ),
            "carbon_lineage_blocker": lineage_target.get("reason", "no-lineage-record"),
            "pubchem_status": pubchem_record.get("status", "not-queried"),
            "pubchem_cid": pubchem_record.get("pubchem", {}).get("CID") or pubchem_record.get("cannabisdb_pubchem_cid"),
            "pubchem_reason": pubchem_record.get("reason"),
            "source": compound.get("source"),
            "source_url": compound.get("source_url"),
        })
    records.sort(key=lambda record: (-record["unresolved_carbon_atoms"], record["label"].lower(), record["cannabisdb_id"]))
    review_queue = [
        {
            "rank": rank,
            "cannabisdb_id": record["cannabisdb_id"],
            "label": record["label"],
            "carbon_atom_count": record["carbon_atom_count"],
            "unresolved_carbon_atoms": record["unresolved_carbon_atoms"],
            "identity_status": record["identity_status"],
            "pubchem_status": record["pubchem_status"],
            "carbon_lineage_status": record["carbon_lineage_status"],
            "blocker": record["carbon_lineage_blocker"],
        }
        for rank, record in enumerate(records, 1)
        if record["unresolved_carbon_atoms"] > 0
    ]
    report = {
        "schema": "cannabis-carbon.named-specialty-inventory.v2",
        "selection": "Source CannabisDB names matching the explicit cannabinoid/cannabis specialty-name regular expression; this is a review inventory, not a chemical-class assertion.",
        "source": str(compounds_path),
        "crosswalk_source": str(crosswalk_path),
        "network_source": str(network_path),
        "lineage_source": str(lineage_path) if lineage_path else None,
        "pubchem_source": str(pubchem_path) if pubchem_path else None,
        "record_count": len(records),
        "carbon_atom_count": sum(record.get("carbon_atom_count") or 0 for record in records),
        "identity_status_counts": {status: sum(record["identity_status"] == status for record in records) for status in ("exact", "candidate", "unresolved")},
        "records_with_reaction_participation": sum(bool(record["reaction_participation"]) for record in records),
        "review_queue_count": len(review_queue),
        "review_queue_unresolved_carbon_atoms": sum(row["unresolved_carbon_atoms"] for row in review_queue),
        "review_queue": review_queue,
        "records": records,
        "claim_boundary": "Name matching identifies records for curation. It does not establish that a compound is endogenous, that a reaction produces it in Cannabis, or that any candidate protein is functional.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {key: report[key] for key in ("record_count", "carbon_atom_count", "identity_status_counts", "records_with_reaction_participation")}
