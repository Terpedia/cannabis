"""Validate the no-silent-omission invariants of published carbon artifacts."""

from __future__ import annotations

import json
from pathlib import Path

from rdkit import Chem


def _carbon_indices(smiles: str | None) -> set[int]:
    molecule = Chem.MolFromSmiles(smiles or "")
    return {atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 6} if molecule else set()


def validate_artifacts(atom_audit_path: Path, mapping_path: Path, balance_path: Path, compounds_path: Path, output: Path) -> dict:
    audit = json.loads(atom_audit_path.read_text())
    mapping = json.loads(mapping_path.read_text())
    balance = json.loads(balance_path.read_text())
    compounds = json.loads(compounds_path.read_text()).get("compounds", [])
    failures: list[dict] = []
    expected_total = 0
    allowed_atom_statuses = {"supported", "candidate", "inferred", "unresolved"}
    observed_status_counts = {status: 0 for status in allowed_atom_statuses}
    audit_records = audit.get("compounds", [])
    audit_ids = [row.get("cannabisdb_id") for row in audit_records]
    for compound_id in sorted({value for value in audit_ids if value is not None}):
        if audit_ids.count(compound_id) > 1:
            failures.append({"kind": "duplicate-compound-audit", "compound_id": compound_id})
    for compound in compounds:
        expected = _carbon_indices(compound.get("smiles"))
        expected_total += len(expected)
        record = next((row for row in audit.get("compounds", []) if row.get("cannabisdb_id") == compound.get("id")), None)
        observed = [index for group in (record or {}).get("groups", []) for index in group.get("atom_indices", [])]
        for group in (record or {}).get("groups", []):
            status = group.get("status")
            if status in observed_status_counts:
                observed_status_counts[status] += len(group.get("atom_indices", []))
        if record is not None and record.get("co2_paths") is not None:
            target_entity_atoms = {
                entity_atom
                for group in record.get("groups", [])
                for entity_atom in group.get("entity_atom_indices", [])
                if entity_atom is not None
            }
            for atom_index, path in record.get("co2_paths", {}).items():
                if path is None or not path:
                    continue
                previous = None
                for step_index, step in enumerate(path):
                    required = ("from_entity_id", "from_atom", "to_entity_id", "to_atom", "status", "provenance")
                    if any(not step.get(field) and step.get(field) != 0 for field in required) or step.get("status") not in {"inferred", "candidate"}:
                        failures.append({"kind": "carbon-path-evidence-missing", "compound_id": compound.get("id"), "atom_index": atom_index, "step": step_index})
                    current_start = (step.get("from_entity_id"), step.get("from_atom"))
                    if previous is None:
                        if step.get("from_entity_id") != "chebi:16526":
                            failures.append({"kind": "carbon-path-does-not-start-at-co2", "compound_id": compound.get("id"), "atom_index": atom_index})
                    elif current_start != previous:
                        failures.append({"kind": "carbon-path-discontinuity", "compound_id": compound.get("id"), "atom_index": atom_index, "step": step_index})
                    previous = (step.get("to_entity_id"), step.get("to_atom"))
                if previous and previous[1] not in target_entity_atoms:
                    failures.append({"kind": "carbon-path-target-mismatch", "compound_id": compound.get("id"), "atom_index": atom_index, "target": previous})
        if record is None:
            failures.append({"kind": "missing-compound-audit", "compound_id": compound.get("id")})
        elif len(observed) != len(set(observed)) or set(observed) != expected:
            failures.append({"kind": "carbon-atom-partition-mismatch", "compound_id": compound.get("id"), "expected": sorted(expected), "observed": sorted(observed)})
        elif any(group.get("status") not in allowed_atom_statuses or not group.get("provenance") or not group.get("reason") for group in record.get("groups", [])):
            failures.append({"kind": "atom-evidence-fields-missing", "compound_id": compound.get("id")})
    for row in mapping.get("reactions", []):
        mapped_rows = row.get("mappings", [])
        product_carbons = row.get("product_carbon_atom_count", 0)
        if len(mapped_rows) != product_carbons:
            failures.append({"kind": "missing-reaction-product-carbon-row", "reaction_id": row.get("reaction_id"), "expected": product_carbons, "observed": len(mapped_rows)})
        for mapped in mapped_rows:
            if mapped.get("status") not in {"inferred", "candidate", "ambiguous", "unresolved"}:
                failures.append({"kind": "invalid-reaction-carbon-status", "reaction_id": row.get("reaction_id"), "product_index": mapped.get("product_index"), "product_atom": mapped.get("product_atom")})
            if mapped.get("status") in {"ambiguous", "unresolved"} and not mapped.get("reason"):
                failures.append({"kind": "reaction-carbon-blocker-missing", "reaction_id": row.get("reaction_id"), "product_index": mapped.get("product_index"), "product_atom": mapped.get("product_atom")})
    imbalanced = [row.get("reaction_id") for row in balance.get("reactions", []) if row.get("status") == "imbalanced"]
    failures.extend({"kind": "imbalanced-reaction", "reaction_id": reaction_id} for reaction_id in imbalanced)
    if audit.get("carbon_atoms_total") != expected_total:
        failures.append({"kind": "global-carbon-total-mismatch", "expected": expected_total, "observed": audit.get("carbon_atoms_total")})
    if sum(observed_status_counts.values()) != expected_total:
        failures.append({"kind": "global-carbon-status-accounting-mismatch", "expected": expected_total, "observed": sum(observed_status_counts.values()), "status_counts": observed_status_counts})
    if audit.get("status_counts") is not None and audit.get("status_counts") != observed_status_counts:
        failures.append({"kind": "reported-carbon-status-count-mismatch", "reported": audit.get("status_counts"), "observed": observed_status_counts})
    result = {"schema": "cannabis-carbon.artifact-validation.v1", "valid": not failures, "checks": {"carbon_atom_partition": not any(f["kind"] in ("missing-compound-audit", "carbon-atom-partition-mismatch", "duplicate-compound-audit") for f in failures), "atom_evidence_fields": not any(f["kind"] == "atom-evidence-fields-missing" for f in failures), "global_carbon_accounting": not any(f["kind"] in ("global-carbon-total-mismatch", "global-carbon-status-accounting-mismatch", "reported-carbon-status-count-mismatch") for f in failures), "carbon_path_integrity": not any(f["kind"].startswith("carbon-path-") for f in failures), "reaction_product_carbon_rows": not any(f["kind"] == "missing-reaction-product-carbon-row" for f in failures), "reaction_mapping_classification": not any(f["kind"] in ("invalid-reaction-carbon-status", "reaction-carbon-blocker-missing") for f in failures), "no_imbalanced_reactions": not imbalanced}, "compounds_checked": len(compounds), "carbon_atoms_expected": expected_total, "carbon_atoms_audited": audit.get("carbon_atoms_total"), "carbon_atom_status_counts": observed_status_counts, "reactions_checked": len(mapping.get("reactions", [])), "imbalanced_reaction_count": len(imbalanced), "failures": failures, "claim_boundary": "This gate validates artifact accounting, evidence-field presence, path continuity, and stoichiometric status. It does not establish enzyme function, isotope tracing, or in-vivo flux."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    return result
