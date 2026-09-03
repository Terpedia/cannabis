"""Conservative RDKit atom-level carbon mapping primitives.

This module maps identical substructures only. It intentionally does not claim
that a maximum-common-substructure match is a biochemical reaction: callers
must attach a reaction record and evidence class separately.
"""

from __future__ import annotations

from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import rdFMCS


@dataclass(frozen=True)
class CarbonMapping:
    reactant_atom: int
    product_atom: int
    status: str
    method: str


def map_conserved_carbons(reactant: Chem.Mol, product: Chem.Mol) -> list[CarbonMapping]:
    """Return conservative carbon correspondences using an MCS.

    Only carbon-to-carbon matches are emitted. Product carbons absent from the
    mapping are therefore explicit unresolved carbon-origin gaps.
    """
    result = rdFMCS.FindMCS(
        [reactant, product], atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder, ringMatchesRingOnly=True,
        completeRingsOnly=True, timeout=30,
    )
    if result.canceled or not result.smartsString:
        return []
    query = Chem.MolFromSmarts(result.smartsString)
    rmatch = reactant.GetSubstructMatch(query)
    pmatch = product.GetSubstructMatch(query)
    return [
        CarbonMapping(r, p, "inferred", "rdkit-mcs-conserved-substructure")
        for r, p in zip(rmatch, pmatch)
        if reactant.GetAtomWithIdx(r).GetAtomicNum() == 6
        and product.GetAtomWithIdx(p).GetAtomicNum() == 6
    ]


def unresolved_product_carbons(product: Chem.Mol, mapping: list[CarbonMapping]) -> list[int]:
    mapped = {item.product_atom for item in mapping}
    return [atom.GetIdx() for atom in product.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in mapped]
