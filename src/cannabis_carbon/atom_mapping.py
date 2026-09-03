"""Conservative atom provenance mapping for Terpedia reaction SMILES."""

from __future__ import annotations

from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import rdFMCS


def _molecules(side: str) -> list[Chem.Mol]:
    molecules = []
    for text in side.split("."):
        mol = Chem.MolFromSmiles(text)
        if mol is not None:
            molecules.append(mol)
    return molecules


def _pair_maps(reactant: Chem.Mol, product: Chem.Mol) -> list[tuple[int, int]]:
    result = rdFMCS.FindMCS(
        [reactant, product], atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder, ringMatchesRingOnly=True,
        completeRingsOnly=True, timeout=5,
    )
    if result.canceled or not result.smartsString:
        return []
    query = Chem.MolFromSmarts(result.smartsString)
    rmatch, pmatch = reactant.GetSubstructMatch(query), product.GetSubstructMatch(query)
    return [(r, p) for r, p in zip(rmatch, pmatch)
            if reactant.GetAtomWithIdx(r).GetAtomicNum() == 6
            and product.GetAtomWithIdx(p).GetAtomicNum() == 6
            # A carbon cannot silently acquire a new heavy-atom position.
            # Bond-order changes (e.g. alcohol -> carbonyl) remain allowed.
            and reactant.GetAtomWithIdx(r).GetDegree() == product.GetAtomWithIdx(p).GetDegree()]


def map_reaction_smiles(reaction_smiles: str) -> dict:
    """Map product carbons to conserved reactant carbons.

    This is structural inference, not isotope tracing. New carbons, ambiguous
    matches, malformed sides, and missing reaction SMILES remain explicit.
    """
    if not reaction_smiles or ">>" not in reaction_smiles:
        return {"status": "unresolved", "reason": "missing-or-invalid-reaction-smiles", "mappings": [], "unresolved_product_carbons": []}
    left, right = reaction_smiles.split(">>", 1)
    reactants, products = _molecules(left), _molecules(right)
    mappings = []
    by_product = defaultdict(list)
    for ri, reactant in enumerate(reactants):
        for pi, product in enumerate(products):
            for ra, pa in _pair_maps(reactant, product):
                by_product[(pi, pa)].append((ri, ra))
    for (pi, pa), choices in by_product.items():
        unique = sorted(set(choices))
        mappings.append({"reactant_index": unique[0][0], "reactant_atom": unique[0][1], "product_index": pi, "product_atom": pa, "status": "inferred" if len(unique) == 1 else "ambiguous", "method": "rdkit-mcs-conserved-substructure", "alternatives": [{"reactant_index": ri, "reactant_atom": ra} for ri, ra in unique]})
    mapped = {(m["product_index"], m["product_atom"]) for m in mappings if m["status"] == "inferred"}
    unresolved = [{"product_index": pi, "product_atom": atom.GetIdx()} for pi, product in enumerate(products) for atom in product.GetAtoms() if atom.GetAtomicNum() == 6 and (pi, atom.GetIdx()) not in mapped]
    status = "inferred" if not unresolved and all(m["status"] == "inferred" for m in mappings) else "unresolved"
    return {"status": status, "mappings": mappings, "unresolved_product_carbons": unresolved, "reactant_count": len(reactants), "product_count": len(products)}
