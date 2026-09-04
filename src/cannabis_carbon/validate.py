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
    for compound in compounds:
        expected = _carbon_indices(compound.get("smiles"))
        expected_total += len(expected)
        record = next((row for row in audit.get("compounds", []) if row.get("cannabisdb_id") == compound.get("id")), None)
        observed = [index for group in (record or {}).get("groups", []) for index in group.get("atom_indices", [])]
        if record is None:
            failures.append({"kind": "missing-compound-audit", "compound_id": compound.get("id")})
        elif len(observed) != len(set(observed)) or set(observed) != expected:
            failures.append({"kind": "carbon-atom-partition-mismatch", "compound_id": compound.get("id"), "expected": sorted(expected), "observed": sorted(observed)})
    for row in mapping.get("reactions", []):
        mapped_rows = row.get("mappings", [])
        product_carbons = row.get("product_carbon_atom_count", 0)
        if len(mapped_rows) != product_carbons:
            failures.append({"kind": "missing-reaction-product-carbon-row", "reaction_id": row.get("reaction_id"), "expected": product_carbons, "observed": len(mapped_rows)})
    imbalanced = [row.get("reaction_id") for row in balance.get("reactions", []) if row.get("status") == "imbalanced"]
    failures.extend({"kind": "imbalanced-reaction", "reaction_id": reaction_id} for reaction_id in imbalanced)
    result = {"schema": "cannabis-carbon.artifact-validation.v1", "valid": not failures, "checks": {"carbon_atom_partition": not any(f["kind"] in ("missing-compound-audit", "carbon-atom-partition-mismatch") for f in failures), "reaction_product_carbon_rows": not any(f["kind"] == "missing-reaction-product-carbon-row" for f in failures), "no_imbalanced_reactions": not imbalanced}, "compounds_checked": len(compounds), "carbon_atoms_expected": expected_total, "carbon_atoms_audited": audit.get("carbon_atoms_total"), "reactions_checked": len(mapping.get("reactions", [])), "imbalanced_reaction_count": len(imbalanced), "failures": failures, "claim_boundary": "This gate validates artifact accounting and stoichiometric status. It does not establish enzyme function, isotope tracing, or in-vivo flux."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    return result
