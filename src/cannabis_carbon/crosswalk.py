from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path

from rdkit import Chem


def build_crosswalk(cannabisdb_sdf: Path, terpedia_network: Path, output: Path) -> dict:
    """Crosswalk CannabisDB compounds to Terpedia ChEBI metabolites by InChIKey."""
    by_key = defaultdict(list)
    by_connectivity = defaultdict(list)
    for index, mol in enumerate(Chem.SDMolSupplier(str(cannabisdb_sdf), removeHs=False, sanitize=True)):
        if mol is None:
            continue
        props = {name: mol.GetProp(name) for name in mol.GetPropNames()}
        key = props.get("INCHI_KEY") or Chem.MolToInchiKey(mol)
        record = {"cannabisdb_id": props.get("DATABASE_ID", f"sdf-record:{index + 1}"), "formula": props.get("FORMULA"), "smiles": props.get("SMILES", Chem.MolToSmiles(mol))}
        by_key[key].append(record)
        if key and "-" in key:
            by_connectivity[key.split("-", 1)[0]].append(record)
    with gzip.open(terpedia_network, "rt", encoding="utf-8") as handle:
        network = json.load(handle)
    matches, unmatched, ambiguous = [], [], []
    metabolites = [e for e in network["entities"] if e.get("type") == "metabolite"]
    terpedia_by_connectivity = {}
    for entity in metabolites:
        attrs = entity.get("attributes", {})
        smiles = attrs.get("canonicalSmiles")
        if not smiles:
            unmatched.append({"terpedia_id": entity["id"], "label": entity.get("label"), "reason": "missing-structure"})
            continue
        mol = Chem.MolFromSmiles(smiles)
        key = Chem.MolToInchiKey(mol) if mol else None
        if key and "-" in key:
            terpedia_by_connectivity.setdefault(key.split("-", 1)[0], []).append(entity)
        candidates = by_key.get(key, []) if key else []
        if len(candidates) == 1:
            matches.append({"terpedia_id": entity["id"], "terpedia_label": entity.get("label"), "cannabisdb": candidates[0], "method": "exact-inchikey"})
        elif len(candidates) > 1:
            ambiguous.append({"terpedia_id": entity["id"], "label": entity.get("label"), "candidates": candidates, "reason": "duplicate-cannabisdb-inchikey"})
        else:
            unmatched.append({"terpedia_id": entity["id"], "label": entity.get("label"), "reason": "no-exact-inchikey-match"})
    exact_ids = {record["cannabisdb"]["cannabisdb_id"] for record in matches}
    candidate_matches = []
    candidate_ambiguous = []
    for connectivity_key, records in by_connectivity.items():
        entities = terpedia_by_connectivity.get(connectivity_key, [])
        for record in records:
            if record["cannabisdb_id"] in exact_ids:
                continue
            if len(entities) == 1:
                candidate_matches.append({"terpedia_id": entities[0]["id"], "terpedia_label": entities[0].get("label"), "cannabisdb": record, "method": "connectivity-inchikey-candidate", "identity_status": "candidate"})
            elif len(entities) > 1:
                candidate_ambiguous.append({"cannabisdb": record, "terpedia_candidates": [{"terpedia_id": e["id"], "label": e.get("label")} for e in entities], "method": "connectivity-inchikey-candidate", "reason": "multiple-terpedia-connectivity-matches"})
    all_cannabisdb_ids = {record["cannabisdb_id"] for records in by_key.values() for record in records}
    matched_cannabisdb_ids = {record["cannabisdb"]["cannabisdb_id"] for record in matches}
    ambiguous_cannabisdb_ids = {candidate["cannabisdb_id"] for record in ambiguous for candidate in record["candidates"]}
    cannabisdb_unmatched_ids = sorted(all_cannabisdb_ids - matched_cannabisdb_ids - ambiguous_cannabisdb_ids)
    report = {"schema": "cannabis-carbon.identity-crosswalk.v1", "method": "RDKit InChIKey derived from exact structures", "cannabisdb_compounds": len(all_cannabisdb_ids), "terpedia_metabolites": len(metabolites), "exact_matches": len(matches), "ambiguous": len(ambiguous), "unmatched": len(unmatched), "terpedia_unmatched": len(unmatched), "cannabisdb_unmatched": len(cannabisdb_unmatched_ids), "connectivity_candidate_matches": len(candidate_matches), "connectivity_candidate_ambiguous": len(candidate_ambiguous), "matches": matches, "candidate_matches": candidate_matches, "candidate_ambiguous_records": candidate_ambiguous, "ambiguous_records": ambiguous, "unmatched_records": unmatched, "cannabisdb_unmatched_ids": cannabisdb_unmatched_ids, "claim_boundary": "Exact matches require full InChIKey identity. Connectivity candidates may differ in stereochemistry, protonation, tautomer, or isotopic state and must not be treated as exact identity or biosynthetic evidence."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {k: report[k] for k in ("cannabisdb_compounds", "terpedia_metabolites", "exact_matches", "ambiguous", "unmatched", "terpedia_unmatched", "cannabisdb_unmatched", "connectivity_candidate_matches", "connectivity_candidate_ambiguous")}
