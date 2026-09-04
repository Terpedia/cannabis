from __future__ import annotations

import json
import re
import unicodedata
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from .terpedia import load_network


_TAUTOMER_ENUMERATOR = rdMolStandardize.TautomerEnumerator()


def _tautomer_key(mol: Chem.Mol | None) -> str | None:
    """Return a candidate-only key after RDKit canonical tautomerization."""
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(_TAUTOMER_ENUMERATOR.Canonicalize(mol))
    except Exception:
        return None


def _name_key(value: str | None) -> str:
    """Normalize a name for candidate lookup without treating it as identity."""
    return re.sub(r"[^a-z0-9]", "", unicodedata.normalize("NFKD", value or "").encode("ascii", "ignore").decode().lower())


def build_crosswalk(cannabisdb_sdf: Path, terpedia_network: Path, output: Path, cannabisdb_catalog: Path | None = None, pubchem_chebi_path: Path | None = None) -> dict:
    """Crosswalk CannabisDB compounds to Terpedia ChEBI metabolites by InChIKey."""
    by_key = defaultdict(list)
    by_connectivity = defaultdict(list)
    by_tautomer = defaultdict(list)
    for index, mol in enumerate(Chem.SDMolSupplier(str(cannabisdb_sdf), removeHs=False, sanitize=True)):
        if mol is None:
            continue
        props = {name: mol.GetProp(name) for name in mol.GetPropNames()}
        key = props.get("INCHI_KEY") or Chem.MolToInchiKey(mol)
        record = {"cannabisdb_id": props.get("DATABASE_ID", f"sdf-record:{index + 1}"), "formula": props.get("FORMULA"), "smiles": props.get("SMILES", Chem.MolToSmiles(mol)), "inchikey": key, "names": sorted({props.get(name, "") for name in ("GENERIC_NAME", "JCHEM_IUPAC", "JCHEM_TRADITIONAL_IUPAC") if props.get(name)})}
        by_key[key].append(record)
        tautomer_key = _tautomer_key(mol)
        if tautomer_key:
            by_tautomer[tautomer_key].append(record)
        if key and "-" in key:
            by_connectivity[key.split("-", 1)[0]].append(record)
    if cannabisdb_catalog:
        catalog = json.loads(cannabisdb_catalog.read_text())
        for row in catalog.get("compounds", []):
            for record in by_key.get(row.get("inchikey"), []):
                record["external_ids"] = row.get("external_ids", {})

    network = load_network(terpedia_network)
    matches, unmatched, ambiguous = [], [], []
    metabolites = [e for e in network["entities"] if e.get("type") == "metabolite"]
    terpedia_by_chebi = defaultdict(list)
    terpedia_by_connectivity = {}
    terpedia_by_tautomer = defaultdict(list)
    terpedia_by_name = defaultdict(list)
    for entity in metabolites:
        attrs = entity.get("attributes", {})
        smiles = attrs.get("canonicalSmiles")
        if not smiles:
            unmatched.append({"terpedia_id": entity["id"], "label": entity.get("label"), "reason": "missing-structure"})
            continue
        mol = Chem.MolFromSmiles(smiles)
        key = Chem.MolToInchiKey(mol) if mol else None
        tautomer_key = _tautomer_key(mol)
        if key and "-" in key:
            terpedia_by_connectivity.setdefault(key.split("-", 1)[0], []).append(entity)
        if tautomer_key:
            terpedia_by_tautomer[tautomer_key].append(entity)
        name_key = _name_key(entity.get("label"))
        if name_key:
            terpedia_by_name[name_key].append(entity)
        chebi_id = entity.get("identifiers", {}).get("chebiId")
        if chebi_id:
            terpedia_by_chebi[str(chebi_id).upper().replace("CHEBI:", "")].append(entity)
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
    connectivity_candidate_ids = {row["cannabisdb"]["cannabisdb_id"] for row in candidate_matches}
    connectivity_ambiguous_ids = {row["cannabisdb"]["cannabisdb_id"] for row in candidate_ambiguous}
    tautomer_candidate_matches = []
    tautomer_candidate_ambiguous = []
    for tautomer_key, records in by_tautomer.items():
        entities = terpedia_by_tautomer.get(tautomer_key, [])
        for record in records:
            compound_id = record["cannabisdb_id"]
            if compound_id in exact_ids or compound_id in connectivity_candidate_ids or compound_id in connectivity_ambiguous_ids:
                continue
            if len(entities) == 1:
                tautomer_candidate_matches.append({"terpedia_id": entities[0]["id"], "terpedia_label": entities[0].get("label"), "cannabisdb": record, "method": "tautomer-canonical-inchikey-candidate", "identity_status": "candidate", "tautomer_key": tautomer_key})
            elif len(entities) > 1:
                tautomer_candidate_ambiguous.append({"cannabisdb": record, "terpedia_candidates": [{"terpedia_id": e["id"], "label": e.get("label")} for e in entities], "method": "tautomer-canonical-inchikey-candidate", "tautomer_key": tautomer_key, "reason": "multiple-terpedia-tautomer-candidates"})
    all_cannabisdb_ids = {record["cannabisdb_id"] for records in by_key.values() for record in records}
    matched_cannabisdb_ids = {record["cannabisdb"]["cannabisdb_id"] for record in matches}
    ambiguous_cannabisdb_ids = {candidate["cannabisdb_id"] for record in ambiguous for candidate in record["candidates"]}
    tautomer_candidate_ids = {row["cannabisdb"]["cannabisdb_id"] for row in tautomer_candidate_matches}
    tautomer_ambiguous_ids = {row["cannabisdb"]["cannabisdb_id"] for row in tautomer_candidate_ambiguous}
    name_candidate_matches = []
    name_candidate_ambiguous = []
    for records in by_key.values():
        for record in records:
            compound_id = record["cannabisdb_id"]
            if compound_id in exact_ids or compound_id in connectivity_candidate_ids or compound_id in connectivity_ambiguous_ids or compound_id in tautomer_candidate_ids or compound_id in tautomer_ambiguous_ids:
                continue
            entities = {entity["id"]: entity for name in record.get("names", []) for entity in terpedia_by_name.get(_name_key(name), [])}
            if len(entities) == 1:
                entity = next(iter(entities.values()))
                terpedia_formula = entity.get("attributes", {}).get("molecularFormula")
                name_candidate_matches.append({"terpedia_id": entity["id"], "terpedia_label": entity.get("label"), "cannabisdb": record, "method": "unique-name-candidate", "identity_status": "candidate", "formula_status": "exact" if not record.get("formula") or not terpedia_formula or terpedia_formula == record["formula"] else "mismatch", "terpedia_formula": terpedia_formula})
            elif len(entities) > 1:
                name_candidate_ambiguous.append({"cannabisdb": record, "terpedia_candidates": [{"terpedia_id": entity["id"], "label": entity.get("label"), "formula": entity.get("attributes", {}).get("molecularFormula")} for entity in entities.values()], "method": "unique-name-candidate", "reason": "multiple-terpedia-name-candidates"})
    name_candidate_ids = {row["cannabisdb"]["cannabisdb_id"] for row in name_candidate_matches}
    name_ambiguous_ids = {row["cannabisdb"]["cannabisdb_id"] for row in name_candidate_ambiguous}
    identifier_matches = []
    identifier_conflicts = []
    if cannabisdb_catalog:
        for records in by_key.values():
            for record in records:
                chebi = record.get("external_ids", {}).get("chebi")
                if not chebi:
                    continue
                entities = terpedia_by_chebi.get(str(chebi).upper().replace("CHEBI:", ""), [])
                if len(entities) != 1:
                    continue
                entity = entities[0]
                terpedia_smiles = entity.get("attributes", {}).get("canonicalSmiles")
                terpedia_mol = Chem.MolFromSmiles(terpedia_smiles) if terpedia_smiles else None
                terpedia_key = Chem.MolToInchiKey(terpedia_mol) if terpedia_mol else None
                row = {"terpedia_id": entity["id"], "terpedia_label": entity.get("label"), "cannabisdb": record, "identifier": f"CHEBI:{chebi}", "method": "cannabisdb-xml-chebi-id"}
                if terpedia_key and terpedia_key == record.get("inchikey"):
                    row["identity_status"] = "structure-verified"
                    identifier_matches.append(row)
                else:
                    row["identity_status"] = "conflict"
                    row["terpedia_inchikey"] = terpedia_key
                    identifier_conflicts.append(row)
    pubchem_chebi_candidate_matches = []
    pubchem_chebi_candidate_ambiguous = []
    pubchem_candidate_ids = set()
    pubchem_ambiguous_ids = set()
    if pubchem_chebi_path and pubchem_chebi_path.exists():
        pubchem_xrefs = json.loads(pubchem_chebi_path.read_text())
        exact_or_candidate_ids = exact_ids
        identifier_record_ids = {row["cannabisdb"]["cannabisdb_id"] for row in identifier_matches + identifier_conflicts}
        for xref in pubchem_xrefs.get("records", []):
            compound_id = xref.get("cannabisdb_id")
            if not compound_id or compound_id in exact_or_candidate_ids or compound_id in identifier_record_ids:
                continue
            entities = {entity["id"]: entity for chebi in xref.get("pubchem_chebi_ids", []) for entity in terpedia_by_chebi.get(str(chebi).upper().replace("CHEBI:", ""), [])}
            if len(entities) == 1:
                entity = next(iter(entities.values()))
                record = next((record for records in by_key.values() for record in records if record.get("cannabisdb_id") == compound_id), None)
                if record:
                    pair = (compound_id, entity["id"])
                    if not any((row["cannabisdb"]["cannabisdb_id"], row["terpedia_id"]) == pair for row in candidate_matches + identifier_matches + identifier_conflicts):
                        pubchem_chebi_candidate_matches.append({"terpedia_id": entity["id"], "terpedia_label": entity.get("label"), "cannabisdb": record, "pubchem_cid": xref.get("cid"), "pubchem_chebi_ids": xref.get("pubchem_chebi_ids", []), "method": "pubchem-chebi-xref-candidate", "identity_status": "candidate", "claim_boundary": "PubChem ChEBI cross-reference is retained as candidate identity evidence; structure, stereochemistry, and reaction context remain unresolved."})
                    pubchem_candidate_ids.add(compound_id)
            elif len(entities) > 1:
                record = next((record for records in by_key.values() for record in records if record.get("cannabisdb_id") == compound_id), None)
                if record:
                    pubchem_chebi_candidate_ambiguous.append({"cannabisdb": record, "pubchem_cid": xref.get("cid"), "pubchem_chebi_ids": xref.get("pubchem_chebi_ids", []), "terpedia_candidates": [{"terpedia_id": entity["id"], "label": entity.get("label")} for entity in entities.values()], "method": "pubchem-chebi-xref-candidate", "reason": "multiple-terpedia-chebi-targets"})
                    pubchem_ambiguous_ids.add(compound_id)
    candidate_matches.extend(pubchem_chebi_candidate_matches)
    candidate_ambiguous.extend(pubchem_chebi_candidate_ambiguous)
    cannabisdb_unmatched_ids = sorted(all_cannabisdb_ids - matched_cannabisdb_ids - ambiguous_cannabisdb_ids - connectivity_candidate_ids - connectivity_ambiguous_ids - tautomer_candidate_ids - tautomer_ambiguous_ids)
    cannabisdb_unmatched_ids = sorted(all_cannabisdb_ids - matched_cannabisdb_ids - ambiguous_cannabisdb_ids - connectivity_candidate_ids - connectivity_ambiguous_ids - tautomer_candidate_ids - tautomer_ambiguous_ids - name_candidate_ids - name_ambiguous_ids - pubchem_candidate_ids - pubchem_ambiguous_ids)
    report = {"schema": "cannabis-carbon.identity-crosswalk.v1", "method": "RDKit InChIKey with candidate-only tautomer, unique-name, CannabisDB XML ChEBI, and PubChem ChEBI cross-reference layers", "cannabisdb_compounds": len(all_cannabisdb_ids), "terpedia_metabolites": len(metabolites), "exact_matches": len(matches), "ambiguous": len(ambiguous), "unmatched": len(unmatched), "terpedia_unmatched": len(unmatched), "cannabisdb_unmatched": len(cannabisdb_unmatched_ids), "connectivity_candidate_matches": len(candidate_matches), "connectivity_candidate_ambiguous": len(candidate_ambiguous), "tautomer_candidate_matches": len(tautomer_candidate_matches), "tautomer_candidate_ambiguous": len(tautomer_candidate_ambiguous), "name_candidate_matches": len(name_candidate_matches), "name_candidate_ambiguous": len(name_candidate_ambiguous), "identifier_matches": len(identifier_matches), "identifier_conflicts": len(identifier_conflicts), "pubchem_chebi_candidate_matches": len(pubchem_chebi_candidate_matches), "pubchem_chebi_candidate_ambiguous": len(pubchem_chebi_candidate_ambiguous), "matches": matches, "identifier_match_records": identifier_matches, "identifier_conflict_records": identifier_conflicts, "candidate_matches": candidate_matches + tautomer_candidate_matches + name_candidate_matches, "candidate_ambiguous_records": candidate_ambiguous + tautomer_candidate_ambiguous + name_candidate_ambiguous, "ambiguous_records": ambiguous, "unmatched_records": unmatched, "cannabisdb_unmatched_ids": cannabisdb_unmatched_ids, "claim_boundary": "Exact matches require full InChIKey identity. XML and PubChem identifier links are reported separately and are structure-verified only when the linked Terpedia structure has the same full InChIKey; conflicts remain unresolved and are not biosynthetic evidence."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {k: report[k] for k in ("cannabisdb_compounds", "terpedia_metabolites", "exact_matches", "ambiguous", "unmatched", "terpedia_unmatched", "cannabisdb_unmatched", "connectivity_candidate_matches", "connectivity_candidate_ambiguous", "tautomer_candidate_matches", "tautomer_candidate_ambiguous", "name_candidate_matches", "name_candidate_ambiguous")}
