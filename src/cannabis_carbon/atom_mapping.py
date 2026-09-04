"""Conservative atom provenance mapping for Terpedia reaction SMILES."""

from __future__ import annotations

from functools import lru_cache

from rdkit import Chem
from rdkit.Chem import rdFMCS
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


@lru_cache(maxsize=4096)
def _mcs_carbon_candidates_cached(reactant_smiles: str, product_smiles: str, product_atom: int) -> tuple[tuple[int, int], ...]:
    """Return unique MCS-derived reactant candidates for one product carbon.

    MCS is used only as a fallback for rearrangements where local topology
    changes. A candidate is retained only when the same MCS atom maps to the
    requested product atom; multiple MCS matches remain ambiguous.
    """
    reactant = Chem.MolFromSmiles(reactant_smiles)
    product = Chem.MolFromSmiles(product_smiles)
    if reactant is None or product is None or not any(atom.GetAtomicNum() == 6 for atom in reactant.GetAtoms()):
        return ()
    result = rdFMCS.FindMCS(
        [reactant, product],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        # Bond order may change in redox and isomerization reactions; element
        # and connectivity conservation remain the constraints here.
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=False,
        completeRingsOnly=False,
        timeout=1,
    )
    if result.canceled or result.numAtoms < 2:
        return ()
    query = Chem.MolFromSmarts(result.smartsString)
    if query is None:
        return ()
    candidates = set()
    for reactant_match in reactant.GetSubstructMatches(query, uniquify=True):
        for product_match in product.GetSubstructMatches(query, uniquify=True):
            if product_atom not in product_match:
                continue
            query_index = product_match.index(product_atom)
            reactant_atom = reactant_match[query_index]
            if reactant.GetAtomWithIdx(reactant_atom).GetAtomicNum() == 6:
                candidates.add((0, reactant_atom))
    return tuple(sorted(candidates))


def _mcs_carbon_candidates(reactants: list[Chem.Mol], product: Chem.Mol, product_atom: int) -> set[tuple[int, int]]:
    candidates = set()
    for ri, reactant in enumerate(reactants):
        pair_candidates = _mcs_carbon_candidates_cached(Chem.MolToSmiles(reactant), Chem.MolToSmiles(product), product_atom)
        candidates.update((ri, atom) for _, atom in pair_candidates)
    return candidates


def _mcs_pair_carbon_candidates(reactants: list[Chem.Mol], products: list[Chem.Mol]) -> dict[tuple[int, int], set[tuple[int, int]]]:
    """Return conservative MCS candidates for every reactant/product pair.

    Local carbon neighborhoods are insufficient for bond-forming and bond-
    breaking reactions.  For multi-reactant reactions we retain pairwise MCS
    correspondences as *candidate* alternatives only: the MCS does not decide
    which substrate supplied a carbon when several correspondences are
    chemically possible.
    """
    candidates = {}
    for pi, product in enumerate(products):
        for atom in product.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            key = (pi, atom.GetIdx())
            pairs = set()
            for ri, reactant in enumerate(reactants):
                for _, ra in _mcs_carbon_candidates_cached(Chem.MolToSmiles(reactant), Chem.MolToSmiles(product), atom.GetIdx()):
                    if reactant.GetAtomWithIdx(ra).GetAtomicNum() == 6:
                        pairs.add((ri, ra))
            if pairs:
                candidates[key] = pairs
    return candidates


def _full_carbon_mcs_mapping(reactants: list[Chem.Mol], products: list[Chem.Mol]) -> list[dict] | None:
    """Map all carbons for a unique one-substrate/one-product MCS."""
    carbon_reactants = [(i, mol) for i, mol in enumerate(reactants) if any(a.GetAtomicNum() == 6 for a in mol.GetAtoms())]
    carbon_products = [(i, mol) for i, mol in enumerate(products) if any(a.GetAtomicNum() == 6 for a in mol.GetAtoms())]
    if len(carbon_reactants) != 1 or len(carbon_products) != 1:
        return None
    ri, reactant = carbon_reactants[0]
    pi, product = carbon_products[0]
    reactant_carbons = [a.GetIdx() for a in reactant.GetAtoms() if a.GetAtomicNum() == 6]
    product_carbons = [a.GetIdx() for a in product.GetAtoms() if a.GetAtomicNum() == 6]
    if len(reactant_carbons) != len(product_carbons):
        return None
    result = rdFMCS.FindMCS([reactant, product], atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareAny, ringMatchesRingOnly=False, completeRingsOnly=False, timeout=2)
    if result.canceled or result.numAtoms < len(product_carbons):
        return None
    query = Chem.MolFromSmarts(result.smartsString)
    if query is None:
        return None
    reactant_matches = reactant.GetSubstructMatches(query, uniquify=True)
    product_matches = product.GetSubstructMatches(query, uniquify=True)
    if len(reactant_matches) != 1 or len(product_matches) != 1:
        return None
    reactant_match, product_match = reactant_matches[0], product_matches[0]
    if set(product_carbons) - set(product_match) or set(reactant_carbons) - set(reactant_match):
        return None
    reactant_by_query = {query_index: atom_index for query_index, atom_index in enumerate(reactant_match)}
    product_by_query = {query_index: atom_index for query_index, atom_index in enumerate(product_match)}
    return [{"product_index": pi, "product_atom": product_atom, "reactant_index": ri, "reactant_atom": reactant_by_query[next(index for index, atom_index in product_by_query.items() if atom_index == product_atom)], "method": "rdkit-full-carbon-mcs-conservation", "status": "inferred"} for product_atom in product_carbons]


def _carbon_skeleton(molecule: Chem.Mol) -> tuple[Chem.Mol, list[int]]:
    """Build a carbon-only graph while retaining its source atom indices."""
    skeleton = Chem.RWMol()
    source_indices = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 6:
            source_indices.append(atom.GetIdx())
            skeleton.AddAtom(Chem.Atom(6))
    source_to_skeleton = {source: index for index, source in enumerate(source_indices)}
    for bond in molecule.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if begin in source_to_skeleton and end in source_to_skeleton:
            skeleton.AddBond(source_to_skeleton[begin], source_to_skeleton[end], Chem.BondType.SINGLE)
    return skeleton.GetMol(), source_indices


def _full_carbon_skeleton_mapping(reactants: list[Chem.Mol], products: list[Chem.Mol]) -> list[dict] | None:
    """Map a unique conserved carbon skeleton when heteroatom chemistry changes."""
    carbon_reactants = [(i, mol) for i, mol in enumerate(reactants) if any(a.GetAtomicNum() == 6 for a in mol.GetAtoms())]
    carbon_products = [(i, mol) for i, mol in enumerate(products) if any(a.GetAtomicNum() == 6 for a in mol.GetAtoms())]
    if len(carbon_reactants) != 1 or len(carbon_products) != 1:
        return None
    ri, reactant = carbon_reactants[0]
    pi, product = carbon_products[0]
    reactant_skeleton, reactant_indices = _carbon_skeleton(reactant)
    product_skeleton, product_indices = _carbon_skeleton(product)
    if len(reactant_indices) != len(product_indices):
        return None
    result = rdFMCS.FindMCS([reactant_skeleton, product_skeleton], atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareAny, ringMatchesRingOnly=False, completeRingsOnly=False, timeout=2)
    if result.canceled or result.numAtoms != len(product_indices):
        return None
    query = Chem.MolFromSmarts(result.smartsString)
    if query is None:
        return None
    reactant_matches = reactant_skeleton.GetSubstructMatches(query, uniquify=True)
    product_matches = product_skeleton.GetSubstructMatches(query, uniquify=True)
    if len(reactant_matches) != 1 or len(product_matches) != 1:
        return None
    reactant_match, product_match = reactant_matches[0], product_matches[0]
    reactant_by_query = {query_index: reactant_indices[atom_index] for query_index, atom_index in enumerate(reactant_match)}
    product_by_query = {query_index: product_indices[atom_index] for query_index, atom_index in enumerate(product_match)}
    return [{"product_index": pi, "product_atom": product_atom, "reactant_index": ri, "reactant_atom": reactant_by_query[next(index for index, atom_index in product_by_query.items() if atom_index == product_atom)], "method": "rdkit-carbon-skeleton-mcs", "status": "inferred"} for product_atom in product_indices]


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
    co2_reactant_indices = [ri for ri, molecule in enumerate(reactants) if Chem.MolToSmiles(molecule) == "O=C=O"]
    # A unique whole-product substructure match is stronger than repeated
    # local carbon signatures for cleavage/decarboxylation-like reactions.
    # It conserves only atoms actually present in the product and does not
    # invent provenance for any eliminated reactant atom.
    if len(reactants) == 1 and len(products) == 1:
        matches = reactants[0].GetSubstructMatches(products[0], uniquify=True)
        if len(matches) == 1:
            mappings = [{"product_index": 0, "product_atom": atom.GetIdx(), "reactant_index": 0, "reactant_atom": matches[0][atom.GetIdx()], "method": "rdkit-unique-product-substructure", "status": "inferred"} for _, atom in product_carbons]
            return {"status": "inferred", "mappings": mappings, "unresolved_product_carbons": [], "product_carbon_atom_count": len(product_carbons), "reactant_count": len(reactants), "product_count": len(products)}
    # If no exact retained-product match exists, accept a complete, unique
    # full-carbon MCS even for simple one-to-one reactions. Equal carbon
    # counts and uniqueness prevent assigning provenance across a carbon gain
    # or an interchangeable structural match.
    full_carbon_mcs = _full_carbon_mcs_mapping(reactants, products)
    if full_carbon_mcs is not None:
        return {"status": "inferred", "mappings": full_carbon_mcs, "unresolved_product_carbons": [], "product_carbon_atom_count": len(product_carbons), "reactant_count": len(reactants), "product_count": len(products)}
    # If only heteroatom connectivity changed, compare the complete carbon
    # skeleton independently.  A unique full skeleton match is still a
    # structural inference; symmetry or carbon rearrangement remains
    # unresolved rather than being assigned by atom order.
    full_carbon_skeleton = _full_carbon_skeleton_mapping(reactants, products)
    if full_carbon_skeleton is not None:
        return {"status": "inferred", "mappings": full_carbon_skeleton, "unresolved_product_carbons": [], "product_carbon_atom_count": len(product_carbons), "reactant_count": len(reactants), "product_count": len(products)}
    # Decarboxylation has one additional product, CO2.  A unique match of the
    # organic product maps its retained carbons; if exactly one reactant
    # carbon remains, that carbon is the released CO2 carbon.
    co2_product_indices = [pi for pi, molecule in enumerate(products) if Chem.MolToSmiles(molecule) == "O=C=O"]
    if len(reactants) == 1 and len(products) == 2 and len(co2_product_indices) == 1:
        organic_indices = [pi for pi in range(len(products)) if pi not in co2_product_indices]
        organic_index = organic_indices[0]
        matches = reactants[0].GetSubstructMatches(products[organic_index], uniquify=True)
        if len(matches) == 1:
            retained = {matches[0][atom.GetIdx()] for atom in products[organic_index].GetAtoms()}
            reactant_carbons = {atom.GetIdx() for atom in reactants[0].GetAtoms() if atom.GetAtomicNum() == 6}
            released = reactant_carbons - {idx for idx in retained if idx in reactant_carbons}
            co2_product = products[co2_product_indices[0]]
            co2_carbons = [atom.GetIdx() for atom in co2_product.GetAtoms() if atom.GetAtomicNum() == 6]
            if len(released) == 1 and len(co2_carbons) == 1:
                mappings = [{"product_index": organic_index, "product_atom": atom.GetIdx(), "reactant_index": 0, "reactant_atom": matches[0][atom.GetIdx()], "method": "rdkit-decarboxylation-substructure", "status": "inferred"} for atom in products[organic_index].GetAtoms() if atom.GetAtomicNum() == 6]
                mappings.append({"product_index": co2_product_indices[0], "product_atom": co2_carbons[0], "reactant_index": 0, "reactant_atom": next(iter(released)), "method": "rdkit-decarboxylation-released-carbon", "status": "inferred"})
                return {"status": "inferred", "mappings": mappings, "unresolved_product_carbons": [], "product_carbon_atom_count": len(product_carbons), "reactant_count": len(reactants), "product_count": len(products)}
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
    # For reactions with multiple substrates, supplement local-signature
    # results with pairwise MCS candidates.  These remain ambiguous/candidate
    # provenance rather than inferred mappings because substrate origin is not
    # uniquely determined by an unlabelled structural MCS.
    # Keep the catalog-wide audit bounded: large cofactors and macromolecular
    # participants make pairwise MCS both expensive and less discriminating.
    # The small-reaction path covers central-carbon rearrangements while the
    # conservative signature mapper remains the fallback for larger records.
    if len(reactants) > 1 and sum(m.GetNumAtoms() for m in reactants + products) <= 60:
        mcs_candidates = _mcs_pair_carbon_candidates(reactants, products)
        for mapping in mappings:
            key = (mapping["product_index"], mapping["product_atom"])
            choices = mcs_candidates.get(key, set())
            if not choices:
                continue
            existing = {(a["reactant_index"], a["reactant_atom"]) for a in mapping.get("alternatives", []) if a.get("reactant_index") is not None}
            if mapping.get("reactant_index") is not None:
                existing.add((mapping["reactant_index"], mapping["reactant_atom"]))
            merged = sorted(existing | choices)
            if mapping["status"] == "inferred":
                continue
            if len(merged) == 1:
                mapping["reactant_index"], mapping["reactant_atom"] = merged[0]
                mapping["status"] = "candidate"
                mapping.pop("reason", None)
                mapping.pop("alternatives", None)
                mapping["method"] = "rdkit-mcs-carbon-candidate"
            else:
                mapping["status"] = "ambiguous"
                mapping["reason"] = "multiple-conserved-reactant-candidates"
                mapping["alternatives"] = [{"reactant_index": ri, "reactant_atom": ra} for ri, ra in merged]
                mapping["method"] = "rdkit-mcs-carbon-candidate"
    # Local signatures intentionally reject many bond-order/ring rearrangements.
    # Use a unique maximum-common-substructure correspondence as a second,
    # still structural, inference method before considering CO2 transfer.
    simple_mcs = len(reactants) == 1 and len(products) == 1 and sum(m.GetNumAtoms() for m in reactants + products) <= 80
    for mapping in mappings:
        if mapping["status"] != "ambiguous" or co2_reactant_indices or not simple_mcs:
            continue
        product = products[mapping["product_index"]]
        candidates = _mcs_carbon_candidates(reactants, product, mapping["product_atom"])
        available = candidates - used_reactant_carbons
        if len(available) == 1:
            ri, ra = next(iter(available))
            used_reactant_carbons.add((ri, ra))
            mapping.update({"reactant_index": ri, "reactant_atom": ra, "status": "inferred", "method": "rdkit-mcs-carbon-conservation-relaxed-bond"})
            mapping.pop("reason", None)
            mapping.pop("alternatives", None)
    # CO2 fixation creates a carbon whose local neighborhood is not conserved.
    # Permit only this explicit inorganic-carbon transfer; do not generalize it
    # to arbitrary unresolved carbons or to CO2 appearing only as a product.
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
