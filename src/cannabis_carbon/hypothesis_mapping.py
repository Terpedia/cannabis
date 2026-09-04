"""RDKit carbon correspondence audit for Terpedia hypothesis connections."""

from __future__ import annotations

import json
from collections import Counter
from rdkit import Chem
from pathlib import Path



def _bounded_pair_mapping(substrate_smiles: str | None, product_smiles: str | None) -> dict:
    """Use only a unique exact product substructure; never run open-ended MCS."""
    substrate = Chem.MolFromSmiles(substrate_smiles or "")
    product = Chem.MolFromSmiles(product_smiles or "")
    if substrate is None or product is None:
        return {"status": "unresolved", "mappings": [], "unresolved_product_carbons": []}
    product_carbons = [atom.GetIdx() for atom in product.GetAtoms() if atom.GetAtomicNum() == 6]
    matches = substrate.GetSubstructMatches(product, uniquify=True)
    if len(matches) != 1:
        return {"status": "candidate" if matches else "unresolved", "mappings": [], "unresolved_product_carbons": product_carbons}
    match = matches[0]
    mappings = [{"product_index": 0, "product_atom": pi, "reactant_index": 0, "reactant_atom": match[pi], "method": "rdkit-unique-product-substructure", "status": "inferred"} for pi in product_carbons]
    return {"status": "inferred", "mappings": mappings, "unresolved_product_carbons": []}


def build_hypothesis_mapping(source: Path, output: Path) -> dict:
    payload = json.loads(source.read_text())
    connections = payload.get("connections", [])

    def map_connection(connection):
        substrate = connection.get("substrate_identity_smiles")
        product = connection.get("product_identity_smiles")
        mapping = _bounded_pair_mapping(substrate, product)
        carbon_mappings = [item for item in mapping.get("mappings", []) if item.get("product_index") == 0]
        status = mapping.get("status", "unresolved")
        return {
            "reaction_id": connection.get("reaction_id"),
            "substrate_terpene_id": connection.get("normalized_substrate_terpene_id"),
            "product_terpene_id": connection.get("normalized_product_terpene_id"),
            "source_type": connection.get("source_type"),
            "evidence_type": connection.get("evidence_type"),
            "substrate_inchikey": connection.get("substrate_inchikey"),
            "product_inchikey": connection.get("product_inchikey"),
            "substrate_carbon_count": int(connection.get("substrate_carbon_count") or 0),
            "product_carbon_count": int(connection.get("product_carbon_count") or 0),
            "mapping_status": status,
            "mapping_method": sorted({item.get("method") for item in carbon_mappings if item.get("method")}),
            "mapped_product_carbon_count": len(carbon_mappings),
            "unresolved_product_carbons": mapping.get("unresolved_product_carbons", []),
            "carbon_mappings": carbon_mappings,
            "claim_boundary": "This is an RDKit structural correspondence for a hypothesis edge, not isotope tracing or proof of in-vivo carbon transfer.",
        }

    rows = [map_connection(connection) for connection in connections]
    report = {
        "schema": "cannabis-carbon.terpedia-hypothesis-carbon-mapping.v1",
        "source": str(source),
        "connection_count": len(rows),
        "summary": dict(sorted(Counter(row["mapping_status"] for row in rows).items())),
        "mapped_product_carbon_atoms": sum(row["mapped_product_carbon_count"] for row in rows),
        "unresolved_product_carbon_atoms": sum(len(row["unresolved_product_carbons"]) for row in rows),
        "reactions": rows,
        "claim_boundary": "This report preserves structural atom correspondences for source-directed Terpedia hypotheses. It does not promote these edges into the balanced reaction network or CO₂ lineage.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {"connections": len(rows), "mapping_status": report["summary"], "mapped_product_carbon_atoms": report["mapped_product_carbon_atoms"], "unresolved_product_carbon_atoms": report["unresolved_product_carbon_atoms"]}
