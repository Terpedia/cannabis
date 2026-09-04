from __future__ import annotations

import json
from pathlib import Path
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from .terpedia import load_network


def _structure_formula(entity: dict) -> tuple[dict[str, int] | None, int | None]:
    smiles = entity.get("attributes", {}).get("canonicalSmiles")
    if not smiles:
        return None, None
    # Rhea generic compounds use wildcard atoms (for example [1*]) and
    # molecular formulas containing R.  Their atom counts are not concrete
    # enough for a stoichiometric balance claim.
    formula_text = entity.get("attributes", {}).get("molecularFormula", "")
    if "*" in smiles or "R" in formula_text:
        return None, entity.get("attributes", {}).get("formalCharge")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    formula = rdMolDescriptors.CalcMolFormula(mol)
    import re
    elements = {symbol: int(count or 1) for symbol, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}
    return elements, Chem.GetFormalCharge(mol)


def _computed_balance(reaction_id: str, entities: dict, statements: list[dict]) -> tuple[dict | None, dict | None]:
    totals = {"left": defaultdict(int), "right": defaultdict(int)}
    charges = {"left": 0, "right": 0}
    complete = True
    participant_count = 0
    for statement in statements:
        if statement.get("subjectId") != reaction_id or statement.get("predicate") not in ("has_reactant", "has_product"):
            continue
        participant_count += 1
        entity = entities.get(statement.get("objectEntityId"), {})
        formula, charge = _structure_formula(entity)
        attrs = entity.get("attributes", {})
        if formula is None:
            formula_text = attrs.get("molecularFormula")
            if formula_text and "R" not in formula_text and all(ch.isalpha() or ch.isdigit() for ch in formula_text):
                import re
                formula = {symbol: int(count or 1) for symbol, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula_text)}
            charge = attrs.get("formalCharge") if attrs.get("formalCharge") is not None else charge
        coefficient = (statement.get("qualifiers") or {}).get("stoichiometricCoefficient", 1)
        side = "left" if statement["predicate"] == "has_reactant" else "right"
        if formula is None or charge is None:
            complete = False
            continue
        for element, count in formula.items(): totals[side][element] += coefficient * count
        charges[side] += coefficient * charge
    if not complete or participant_count == 0 or not totals["left"] or not totals["right"]:
        return None, None
    differences = {element: totals["right"][element] - totals["left"][element] for element in sorted(set(totals["left"]) | set(totals["right"])) if totals["right"][element] != totals["left"][element]}
    return {"status": "balanced" if not differences else "imbalanced", "differences": differences}, {"status": "balanced" if charges["left"] == charges["right"] else "imbalanced", "left": charges["left"], "right": charges["right"], "difference": charges["right"] - charges["left"]}


def audit_balances(network_path: Path, output: Path) -> dict:
    network = load_network(network_path)
    rows = []
    entities = {e["id"]: e for e in network["entities"]}
    for entity in (e for e in network["entities"] if e.get("type") == "biochemical_reaction"):
        attrs = entity.get("attributes", {})
        element = attrs.get("elementBalance") or {"status": "missing"}
        charge = attrs.get("chargeBalance") or {"status": "missing"}
        computed_element, computed_charge = _computed_balance(entity["id"], entities, network["statements"])
        computed_status = "balanced" if computed_element and computed_charge and computed_element["status"] == "balanced" and computed_charge["status"] == "balanced" else "imbalanced" if computed_element and computed_charge else "not_auditable"
        rows.append({"reaction_id": entity["id"], "label": entity.get("label"), "equation": attrs.get("equation"), "element_balance": element, "charge_balance": charge, "computed_element_balance": computed_element, "computed_charge_balance": computed_charge, "status": "balanced" if element.get("status") == "balanced" and charge.get("status") == "balanced" else computed_status})
    report = {"schema": "cannabis-carbon.phase1-balance-audit.v1", "phase": "Phase 1", "reaction_count": len(rows), "summary": {"fully_balanced": sum(r["status"] == "balanced" for r in rows), "imbalanced": sum(r["status"] == "imbalanced" for r in rows), "not_auditable": sum(r["status"] == "not_auditable" for r in rows), "element_balanced": sum(r["element_balance"].get("status") == "balanced" for r in rows), "charge_balanced": sum(r["charge_balance"].get("status") == "balanced" for r in rows), "computed_fully_balanced": sum(bool(r["computed_element_balance"] and r["computed_charge_balance"] and r["computed_element_balance"].get("status") == "balanced" and r["computed_charge_balance"].get("status") == "balanced") for r in rows), "computed_imbalanced": sum(bool(r["computed_element_balance"] and r["computed_charge_balance"] and (r["computed_element_balance"].get("status") == "imbalanced" or r["computed_charge_balance"].get("status") == "imbalanced")) for r in rows)}, "reactions": rows, "claim_boundary": "Balance is a Phase 1 stoichiometric validity gate; it does not establish enzyme identity, direction, carbon atom mapping, or in-vivo flux."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return report["summary"]
