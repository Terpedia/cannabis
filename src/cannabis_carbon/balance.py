from __future__ import annotations

import gzip
import json
from pathlib import Path


def audit_balances(network_path: Path, output: Path) -> dict:
    with gzip.open(network_path, "rt", encoding="utf-8") as handle:
        network = json.load(handle)
    rows = []
    for entity in (e for e in network["entities"] if e.get("type") == "biochemical_reaction"):
        attrs = entity.get("attributes", {})
        element = attrs.get("elementBalance") or {"status": "missing"}
        charge = attrs.get("chargeBalance") or {"status": "missing"}
        rows.append({"reaction_id": entity["id"], "label": entity.get("label"), "equation": attrs.get("equation"), "element_balance": element, "charge_balance": charge, "status": "balanced" if element.get("status") == "balanced" and charge.get("status") == "balanced" else "not_auditable" if "not_auditable" in (element.get("status"), charge.get("status")) else "imbalanced"})
    report = {"schema": "cannabis-carbon.phase1-balance-audit.v1", "phase": "Phase 1", "reaction_count": len(rows), "summary": {"fully_balanced": sum(r["status"] == "balanced" for r in rows), "imbalanced": sum(r["status"] == "imbalanced" for r in rows), "not_auditable": sum(r["status"] == "not_auditable" for r in rows), "element_balanced": sum(r["element_balance"].get("status") == "balanced" for r in rows), "charge_balanced": sum(r["charge_balance"].get("status") == "balanced" for r in rows)}, "reactions": rows, "claim_boundary": "Balance is a Phase 1 stoichiometric validity gate; it does not establish enzyme identity, direction, carbon atom mapping, or in-vivo flux."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return report["summary"]
