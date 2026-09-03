"""Conservative atom provenance mapping for Terpedia reaction SMILES."""

from __future__ import annotations

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.warning")


def _molecules(side: str) -> list[Chem.Mol]:
    molecules = []
    for text in side.split("."):
        mol = Chem.MolFromSmiles(text)
        if mol is not None:
            molecules.append(mol)
    return molecules


def _carbon_signature(atom: Chem.Atom) -> tuple:
    """A bond-order-independent local invariant for conservative mapping."""
    return (atom.GetDegree(), atom.IsInRing(), atom.GetIsAromatic(), tuple(sorted(n.GetAtomicNum() for n in atom.GetNeighbors())))


def map_reaction_smiles(reaction_smiles: str) -> dict:
    """Map product carbons to conserved reactant carbons.

    This is structural inference, not isotope tracing. New carbons, ambiguous
    matches, malformed sides, and missing reaction SMILES remain explicit.
    """
    if not reaction_smiles or ">>" not in reaction_smiles:
        return {"status": "unresolved", "reason": "missing-or-invalid-reaction-smiles", "mappings": [], "unresolved_product_carbons": [], "product_carbon_atom_count": 0}
    left, right = reaction_smiles.split(">>", 1)
    reactants, products = _molecules(left), _molecules(right)
    mappings = []
    reactant_carbons = [(ri, atom.GetIdx(), _carbon_signature(atom)) for ri, mol in enumerate(reactants) for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    used_reactant_carbons = set()
    product_carbons = [(pi, atom) for pi, product in enumerate(products) for atom in product.GetAtoms() if atom.GetAtomicNum() == 6]
    # Process constrained product atoms first. A reactant atom may be used only
    # once; this prevents repeated local signatures from fabricating carbon
    # conservation. Non-unique matches are retained as explicit ambiguity.
    candidate_rows = []
    for pi, atom in product_carbons:
        choices = [(ri, ra) for ri, ra, signature in reactant_carbons if signature == _carbon_signature(atom)]
        candidate_rows.append((len(choices), pi, atom, choices))
    for _, pi, atom, choices in sorted(candidate_rows, key=lambda row: (row[0], row[1], row[2].GetIdx())):
        base = {"product_index": pi, "product_atom": atom.GetIdx(), "method": "rdkit-carbon-neighborhood"}
        if len(choices) == 1 and choices[0] not in used_reactant_carbons:
            used_reactant_carbons.add(choices[0])
            mappings.append({**base, "reactant_index": choices[0][0], "reactant_atom": choices[0][1], "status": "inferred"})
        elif len(choices) > 1:
            mappings.append({**base, "status": "ambiguous", "reason": "multiple-conserved-reactant-candidates", "alternatives": [{"reactant_index": ri, "reactant_atom": ra} for ri, ra in choices]})
        else:
            reason = "reactant-carbon-already-assigned" if choices else "no-conserved-reactant-carbon"
            mappings.append({**base, "status": "unresolved", "reason": reason})
    # CO2 fixation creates a carbon whose local neighborhood is not conserved.
    # Permit only this explicit inorganic-carbon transfer; do not generalize it
    # to arbitrary unresolved carbons or to CO2 appearing only as a product.
    co2_reactant_indices = [ri for ri, molecule in enumerate(reactants) if Chem.MolToSmiles(molecule) == "O=C=O"]
    for ri in co2_reactant_indices:
        for mapping in sorted((m for m in mappings if m["status"] != "inferred"), key=lambda m: (-sum(n.GetAtomicNum() == 8 for n in products[m["product_index"]].GetAtomWithIdx(m["product_atom"]).GetNeighbors()), m["product_index"], m["product_atom"])):
            product_atom = products[mapping["product_index"]].GetAtomWithIdx(mapping["product_atom"])
            oxygen_neighbors = sum(n.GetAtomicNum() == 8 for n in product_atom.GetNeighbors())
            if oxygen_neighbors < 2:
                continue
            mapping.update({"reactant_index": ri, "reactant_atom": next(atom.GetIdx() for atom in reactants[ri].GetAtoms() if atom.GetAtomicNum() == 6), "status": "inferred", "method": "rdkit-co2-carbon-source"})
            mapping.pop("reason", None)
            mapping.pop("alternatives", None)
            break
        else:
            continue
        break
    # Some fixation reactions produce a carbonyl rather than a carboxylate.
    # Keep this weaker assignment visibly candidate-supported.
    if co2_reactant_indices and any(m["status"] != "inferred" for m in mappings):
        ri = co2_reactant_indices[0]
        for mapping in sorted((m for m in mappings if m["status"] != "inferred"), key=lambda m: (m["product_index"], m["product_atom"])):
            mapping.update({"reactant_index": ri, "reactant_atom": next(atom.GetIdx() for atom in reactants[ri].GetAtoms() if atom.GetAtomicNum() == 6), "status": "candidate", "method": "rdkit-co2-carbon-source-candidate"})
            mapping.pop("reason", None)
            mapping.pop("alternatives", None)
            break
    unresolved = [{"product_index": m["product_index"], "product_atom": m["product_atom"], "status": m["status"], "reason": m.get("reason")} for m in mappings if m["status"] not in ("inferred", "candidate")]
    status = "inferred" if not unresolved and all(m["status"] == "inferred" for m in mappings) else "unresolved"
    return {"status": status, "mappings": mappings, "unresolved_product_carbons": unresolved, "product_carbon_atom_count": sum(atom.GetAtomicNum() == 6 for product in products for atom in product.GetAtoms()), "reactant_count": len(reactants), "product_count": len(products)}
