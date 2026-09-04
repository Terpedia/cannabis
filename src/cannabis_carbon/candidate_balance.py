"""Stoichiometric audit for the Terpedia identity-set candidate expansion."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from .balance import _reaction_smiles_balance


def audit_candidate_expansion_balances(source: Path, output: Path) -> dict:
    payload = json.loads(source.read_text())
    rows = []
    seen = set()
    edge_statuses = Counter()
    for edge in payload.get("rows", []):
        reaction_id = edge.get("reaction_id")
        reaction_smarts = edge.get("reaction_smarts")
        key = (reaction_id, reaction_smarts)
        element, charge = _reaction_smiles_balance(reaction_smarts)
        edge_status = (
            "balanced" if element and charge and element["status"] == "balanced" and charge["status"] == "balanced"
            else "imbalanced" if element and charge
            else "not_auditable"
        )
        edge_statuses[edge_status] += 1
        if key in seen:
            continue
        seen.add(key)
        rows.append({
            "reaction_id": reaction_id,
            "reaction_smarts": reaction_smarts,
            "source_types": sorted({r.get("source_type") for r in payload.get("rows", []) if (r.get("reaction_id"), r.get("reaction_smarts")) == key and r.get("source_type")}),
            "edge_count": sum((r.get("reaction_id"), r.get("reaction_smarts")) == key for r in payload.get("rows", [])),
            "element_balance": element,
            "charge_balance": charge,
            "status": edge_status,
            "claim_boundary": "Balanced status is a stoichiometric screen of a source reaction SMARTS; it does not establish enzyme identity, physiological direction, atom mapping, or in-vivo Cannabis biosynthesis.",
        })
    rows.sort(key=lambda r: (r.get("reaction_id") or "", r.get("reaction_smarts") or ""))
    report = {
        "schema": "cannabis-carbon.identity-set-candidate-expansion-balance-audit.v1",
        "source": str(source),
        "edge_count": len(payload.get("rows", [])),
        "unique_reaction_count": len(rows),
        "edge_status_counts": dict(sorted(edge_statuses.items())),
        "reaction_status_counts": dict(sorted(Counter(row["status"] for row in rows).items())),
        "reactions": rows,
        "claim_boundary": "This is a Phase 1 stoichiometric audit of the separate candidate expansion layer. Only balanced rows are eligible for downstream candidate-path review; no row is promoted to confirmed Cannabis metabolism by this audit.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {key: report[key] for key in ("edge_count", "unique_reaction_count", "edge_status_counts", "reaction_status_counts")}
