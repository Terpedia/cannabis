from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path

from rdkit import Chem


def build_crosswalk(cannabisdb_sdf: Path, terpedia_network: Path, output: Path) -> dict:
    """Crosswalk CannabisDB compounds to Terpedia ChEBI metabolites by InChIKey."""
    by_key = defaultdict(list)
    for index, mol in enumerate(Chem.SDMolSupplier(str(cannabisdb_sdf), removeHs=False, sanitize=True)):
        if mol is None:
            continue
        props = {name: mol.GetProp(name) for name in mol.GetPropNames()}
        key = props.get("INCHI_KEY") or Chem.MolToInchiKey(mol)
        by_key[key].append({"cannabisdb_id": props.get("DATABASE_ID", f"sdf-record:{index + 1}"), "formula": props.get("FORMULA"), "smiles": props.get("SMILES", Chem.MolToSmiles(mol))})
    with gzip.open(terpedia_network, "rt", encoding="utf-8") as handle:
        network = json.load(handle)
    matches, unmatched, ambiguous = [], [], []
    for entity in (e for e in network["entities"] if e.get("type") == "metabolite"):
        attrs = entity.get("attributes", {})
        smiles = attrs.get("canonicalSmiles")
        if not smiles:
            unmatched.append({"terpedia_id": entity["id"], "label": entity.get("label"), "reason": "missing-structure"})
            continue
        mol = Chem.MolFromSmiles(smiles)
        key = Chem.MolToInchiKey(mol) if mol else None
        candidates = by_key.get(key, []) if key else []
        if len(candidates) == 1:
            matches.append({"terpedia_id": entity["id"], "terpedia_label": entity.get("label"), "cannabisdb": candidates[0], "method": "exact-inchikey"})
        elif len(candidates) > 1:
            ambiguous.append({"terpedia_id": entity["id"], "label": entity.get("label"), "candidates": candidates, "reason": "duplicate-cannabisdb-inchikey"})
        else:
            unmatched.append({"terpedia_id": entity["id"], "label": entity.get("label"), "reason": "no-exact-inchikey-match"})
    report = {"schema": "cannabis-carbon.identity-crosswalk.v1", "method": "RDKit InChIKey derived from exact structures", "cannabisdb_compounds": sum(len(v) for v in by_key.values()), "terpedia_metabolites": len([e for e in network["entities"] if e.get("type") == "metabolite"]), "exact_matches": len(matches), "ambiguous": len(ambiguous), "unmatched": len(unmatched), "matches": matches, "ambiguous_records": ambiguous, "unmatched_records": unmatched, "claim_boundary": "An exact structure identity does not establish that a CannabisDB compound is biosynthesized by Cannabis."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {k: report[k] for k in ("cannabisdb_compounds", "terpedia_metabolites", "exact_matches", "ambiguous", "unmatched")}
