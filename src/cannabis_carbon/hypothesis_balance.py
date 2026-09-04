"""Stoichiometric audit for the separate Terpedia hypothesis-edge layer."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from .balance import _reaction_smiles_balance


def audit_hypothesis_balances(source: Path, output: Path) -> dict:
    payload = json.loads(source.read_text())
    rows = []
    for connection in payload.get("connections", []):
        element, charge = _reaction_smiles_balance(connection.get("reaction_smarts"))
        status = (
            "balanced" if element and charge and element["status"] == "balanced" and charge["status"] == "balanced"
            else "imbalanced" if element and charge
            else "not_auditable"
        )
        rows.append({
            "reaction_id": connection.get("reaction_id"),
            "substrate_terpene_id": connection.get("normalized_substrate_terpene_id"),
            "product_terpene_id": connection.get("normalized_product_terpene_id"),
            "source_type": connection.get("source_type"),
            "evidence_type": connection.get("evidence_type"),
            "reaction_smarts": connection.get("reaction_smarts"),
            "element_balance": element,
            "charge_balance": charge,
            "status": status,
            "claim_boundary": "This is a stoichiometric screen of a hypothesis edge; balanced status does not establish reaction direction, enzyme function, or in-vivo Cannabis biosynthesis.",
        })
    report = {
        "schema": "cannabis-carbon.terpedia-hypothesis-balance-audit.v1",
        "source": str(source),
        "connection_count": len(rows),
        "summary": dict(sorted(Counter(row["status"] for row in rows).items())),
        "source_type_summary": {
            source_type: dict(sorted(Counter(row["status"] for row in rows if row["source_type"] == source_type).items()))
            for source_type in sorted({row["source_type"] for row in rows})
        },
        "reactions": rows,
        "claim_boundary": "This audit classifies source-directed hypothesis reaction SMARTS for Phase 1 review. It does not promote any edge into the balanced reaction network or CO₂ lineage without independent validation.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return report["summary"]
