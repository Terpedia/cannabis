"""Extract a focused upstream reaction neighborhood from Terpedia BigQuery."""

from __future__ import annotations

import json
import subprocess
from collections import Counter
from datetime import date
from pathlib import Path

from rdkit import Chem

from .identity_set import SOURCE_TABLE
from .atom_mapping import map_reaction_smiles, map_identity_pair_smiles


EDGE_TABLE = "terpedia-489015.terpedia_core.terpene_metabolic_map_edges_current"


def _carbon_bearing_required_substrates(value: str | None) -> list[dict]:
    """Return carbon-bearing required substrates as unresolved alternatives."""
    try:
        structures = json.loads(value or "[]")
    except json.JSONDecodeError:
        return []
    candidates = []
    for index, structure in enumerate(structures):
        molecule = Chem.MolFromSmiles(structure) if isinstance(structure, str) else None
        if molecule is None:
            continue
        carbon_count = sum(atom.GetAtomicNum() == 6 for atom in molecule.GetAtoms())
        if not carbon_count:
            continue
        candidates.append({
            "input_index": index,
            "smiles": structure,
            "canonical_smiles": Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True),
            "carbon_count": carbon_count,
            "source": "Terpedia required_substrate_structures_json",
            "identity_status": "unresolved-candidate-structure",
        })
    return candidates


def map_identity_set_upstream(input_path: Path, output: Path) -> dict:
    """Add conservative RDKit carbon mappings to a focused upstream queue."""
    report = json.loads(input_path.read_text())
    rows = []
    status_counts = Counter()
    pair_status_counts = Counter()
    pair_mapped_atoms = 0
    pair_unresolved_atoms = 0
    missing_precursor_rows_with_candidates = 0
    candidate_pair_status_counts = Counter()
    candidate_pair_mapped_atoms = 0
    candidate_pair_unresolved_atoms = 0
    for row in report.get("rows", []):
        mapped = map_reaction_smiles(row.get("reaction_smarts"))
        product = report.get("identity_records", {}).get(row.get("product_identity_id"), {}).get("smiles")
        precursor = report.get("identity_records", {}).get(row.get("precursor_identity_id"), {}).get("smiles")
        pair_mapping = map_identity_pair_smiles(precursor, product)
        precursor_candidates = [] if precursor else _carbon_bearing_required_substrates(row.get("required_substrate_structures_json"))
        if precursor_candidates:
            missing_precursor_rows_with_candidates += 1
            for candidate in precursor_candidates:
                candidate_mapping = map_identity_pair_smiles(candidate["canonical_smiles"], product)
                candidate["target_pair_carbon_mapping"] = candidate_mapping
                candidate_pair_status_counts[candidate_mapping.get("status", "unresolved")] += 1
                candidate_pair_mapped_atoms += len(candidate_mapping.get("mappings", []))
                candidate_pair_unresolved_atoms += len(candidate_mapping.get("unresolved_product_carbons", []))
        pair_status_counts[pair_mapping.get("status", "unresolved")] += 1
        pair_mapped_atoms += len(pair_mapping.get("mappings", []))
        pair_unresolved_atoms += len(pair_mapping.get("unresolved_product_carbons", []))
        enriched = {**row, "carbon_mapping": mapped, "target_pair_carbon_mapping": pair_mapping, "precursor_structure_candidates": precursor_candidates}
        rows.append(enriched)
        status_counts[mapped.get("status", "unresolved")] += 1
    mapped_report = {
        **report,
        "schema": "cannabis-carbon.terpedia-identity-set-upstream-mapped.v1",
        "mapping_method": "RDKit map_reaction_smiles applied to each source reaction SMARTS; mappings are structural provenance candidates, not isotope tracing.",
        "mapping_status_counts": dict(sorted(status_counts.items())),
        "target_pair_mapping_status_counts": dict(sorted(pair_status_counts.items())),
        "mapped_product_carbon_atoms": sum(len(row["carbon_mapping"].get("mappings", [])) for row in rows),
        "unresolved_product_carbon_atoms": sum(len(row["carbon_mapping"].get("unresolved_product_carbons", [])) for row in rows),
        "target_pair_mapped_carbon_atoms": pair_mapped_atoms,
        "target_pair_unresolved_carbon_atoms": pair_unresolved_atoms,
        "missing_precursor_rows_with_structure_candidates": missing_precursor_rows_with_candidates,
        "candidate_pair_mapping_status_counts": dict(sorted(candidate_pair_status_counts.items())),
        "candidate_pair_mapped_carbon_atoms": candidate_pair_mapped_atoms,
        "candidate_pair_unresolved_carbon_atoms": candidate_pair_unresolved_atoms,
        "rows": rows,
        "claim_boundary": "RDKit mappings conserve structurally matched carbon atoms and retain unresolved product atoms explicitly. They do not establish reaction direction, enzyme function, isotope tracing, or endogenous Cannabis biosynthesis.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(mapped_report, separators=(",", ":")) + "\n")
    return {key: mapped_report[key] for key in ("edge_count", "mapped_product_carbon_atoms", "unresolved_product_carbon_atoms", "mapping_status_counts")}


def refresh_identity_set_upstream(compounds_path: Path, output: Path, bq_binary: str = "bq") -> dict:
    """Extract current Terpedia edges whose products match CannabisDB identity-set records.

    The result is intentionally a focused evidence queue, not an assertion that
    every extracted edge is a Cannabis biosynthetic reaction.
    """
    catalog = json.loads(compounds_path.read_text())
    compounds = catalog.get("compounds", catalog if isinstance(catalog, list) else [])
    by_key = {c["inchikey"]: c for c in compounds if c.get("inchikey")}
    keys = sorted(by_key)
    literals = ",".join("'" + key.replace("'", "\\'") + "'" for key in keys)
    query = f"""
WITH target AS (
  SELECT DISTINCT terpene_id, inchikey, inchi, smiles, molecular_formula,
         carbon_count, source_crossrefs
  FROM `{SOURCE_TABLE}`
  WHERE inchikey IN ({literals})
),
identity AS (
  SELECT DISTINCT terpene_id, inchikey, inchi, smiles, molecular_formula,
         carbon_count, source_crossrefs
  FROM `{SOURCE_TABLE}`
)
SELECT
  e.product_terpene_id, e.precursor_terpene_id, e.reaction_id,
  e.source_type, e.evidence_type, e.structure_match_mode,
  e.reaction_smarts, e.required_substrate_structures_json,
  e.missing_corpus_substrates_json, e.source_dataset, e.source_url,
  e.source_uniprot_id, e.source_genbank_id, e.source_ec_number,
  e.claim_boundary,
  target.inchikey AS product_inchikey,
  target.inchi AS product_inchi,
  target.smiles AS product_smiles,
  target.molecular_formula AS product_molecular_formula,
  target.carbon_count AS product_carbon_count,
  target.source_crossrefs AS product_source_crossrefs,
  product_identity.inchikey AS precursor_inchikey,
  product_identity.inchi AS precursor_inchi,
  product_identity.smiles AS precursor_smiles,
  product_identity.molecular_formula AS precursor_molecular_formula,
  product_identity.carbon_count AS precursor_carbon_count,
  product_identity.source_crossrefs AS precursor_source_crossrefs
FROM `{EDGE_TABLE}` e
JOIN target ON target.terpene_id = e.product_terpene_id
LEFT JOIN identity AS product_identity
  ON product_identity.terpene_id = e.precursor_terpene_id
"""
    completed = subprocess.run(
        [bq_binary, "query", "--use_legacy_sql=false", "--format=prettyjson", "--max_rows=10000", query],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = json.loads(completed.stdout)
    rows.sort(key=lambda row: (row.get("product_terpene_id") or "", row.get("reaction_id") or "", row.get("precursor_terpene_id") or ""))
    identity_fields = ("inchikey", "inchi", "smiles", "molecular_formula", "carbon_count", "source_crossrefs")
    identities = {}
    product_cannabisdb_ids = {}
    compact_rows = []
    for row in rows:
        product_id = row.get("product_terpene_id")
        precursor_id = row.get("precursor_terpene_id")
        if product_id and product_id not in identities:
            identities[product_id] = {key: row.get(f"product_{key}") for key in identity_fields}
        if precursor_id and precursor_id not in identities and row.get("precursor_inchikey"):
            identities[precursor_id] = {key: row.get(f"precursor_{key}") for key in identity_fields}
        if product_id and row.get("product_inchikey") in by_key:
            product_cannabisdb_ids.setdefault(product_id, []).append(by_key[row["product_inchikey"]]["id"])
        expanded_identity_fields = {f"product_{key}" for key in identity_fields} | {f"precursor_{key}" for key in identity_fields}
        compact = {key: value for key, value in row.items() if key not in expanded_identity_fields}
        compact.update({"product_identity_id": product_id, "precursor_identity_id": precursor_id})
        compact_rows.append(compact)
    for identity_id in identities:
        product_cannabisdb_ids[identity_id] = sorted(set(product_cannabisdb_ids.get(identity_id, [])))
    target_by_terpene = {}
    for row in compact_rows:
        target_by_terpene.setdefault(row.get("product_terpene_id"), [])
        if row.get("product_identity_id") in product_cannabisdb_ids:
            target_by_terpene[row["product_terpene_id"]].extend(product_cannabisdb_ids[row["product_identity_id"]])
    product_keys = {key for key in product_cannabisdb_ids if product_cannabisdb_ids[key]}
    report = {
        "schema": "cannabis-carbon.terpedia-identity-set-upstream.v1",
        "source_table": EDGE_TABLE,
        "identity_set_table": SOURCE_TABLE,
        "retrieved_at": date.today().isoformat(),
        "method": "Current Terpedia upstream edges joined to exact full-InChIKey identity-set matches for the CannabisDB catalog.",
        "cannabisdb_compounds": len(compounds),
        "identity_set_target_count": len(product_keys),
        "edge_count": len(rows),
        "precursor_identity_count": len({row.get("precursor_terpene_id") for row in rows if row.get("precursor_terpene_id")}),
        "reaction_count": len({row.get("reaction_id") for row in rows if row.get("reaction_id")}),
        "source_type_counts": dict(sorted(Counter(row.get("source_type") or "unknown" for row in rows).items())),
        "evidence_type_counts": dict(sorted(Counter(row.get("evidence_type") or "unknown" for row in rows).items())),
        "target_terpene_ids": sorted(target_by_terpene),
        "identity_records": identities,
        "product_cannabisdb_ids": {key: sorted(set(value)) for key, value in product_cannabisdb_ids.items() if value},
        "rows": compact_rows,
        "claim_boundary": "These are source-linked Terpedia reaction edges entering identity-set-linked CannabisDB structures. They are an upstream evidence queue, not proof of endogenous Cannabis biosynthesis, enzyme function, reaction direction, or CO2 carbon provenance.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {key: report[key] for key in ("cannabisdb_compounds", "identity_set_target_count", "edge_count", "precursor_identity_count", "reaction_count")}
