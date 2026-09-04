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
from .terpedia import load_network


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


def build_identity_set_core_bridges(mapped_path: Path, network_path: Path, output: Path) -> dict:
    """Find candidate bridges from core Terpedia metabolites to CDB targets."""
    report = json.loads(mapped_path.read_text())
    network = load_network(network_path)
    core_by_key = {}
    for entity in network.get("entities", []):
        if entity.get("type") != "metabolite":
            continue
        smiles = entity.get("attributes", {}).get("canonicalSmiles")
        molecule = Chem.MolFromSmiles(smiles) if smiles else None
        if molecule is not None:
            core_by_key.setdefault(Chem.MolToInchiKey(molecule), []).append(entity["id"])
    identity_keys = {}
    for identity_id, record in report.get("identity_records", {}).items():
        molecule = Chem.MolFromSmiles(record.get("smiles")) if record.get("smiles") else None
        if molecule is not None:
            identity_keys[identity_id] = Chem.MolToInchiKey(molecule)
    cdb_by_terpene = report.get("product_cannabisdb_ids", {})
    bridges = []
    for row in report.get("rows", []):
        precursor_id = row.get("precursor_identity_id")
        product_id = row.get("product_identity_id")
        core_entities = core_by_key.get(identity_keys.get(precursor_id), [])
        cdb_ids = cdb_by_terpene.get(product_id, [])
        if not core_entities or not cdb_ids:
            continue
        for core_id in core_entities:
            for cdb_id in cdb_ids:
                bridges.append({
                    "cannabisdb_id": cdb_id,
                    "core_entity_id": core_id,
                    "precursor_identity_id": precursor_id,
                    "product_identity_id": product_id,
                    "reaction_id": row.get("reaction_id"),
                    "source_type": row.get("source_type"),
                    "evidence_type": row.get("evidence_type"),
                    "source_url": row.get("source_url"),
                    "reaction_smarts": row.get("reaction_smarts"),
                    "target_pair_carbon_mapping": row.get("target_pair_carbon_mapping"),
                    "status": "candidate",
                    "claim_boundary": "This exact-structure bridge connects a Terpedia core metabolite to a CannabisDB identity-set target through a source-linked reaction edge. It is a candidate hypothesis, not proof of reaction direction, enzyme activity, endogenous biosynthesis, or CO2 carbon provenance.",
                })
    bridges.sort(key=lambda bridge: (bridge["cannabisdb_id"], bridge["reaction_id"] or "", bridge["core_entity_id"]))
    result = {
        "schema": "cannabis-carbon.identity-set-core-bridges.v1",
        "source_mapped_upstream": str(mapped_path),
        "source_network": str(network_path),
        "method": "Exact InChIKey matching of identity-set precursor structures to Terpedia core metabolite structures, retaining source edge and target-pair RDKit evidence.",
        "bridge_count": len(bridges),
        "target_count": len({bridge["cannabisdb_id"] for bridge in bridges}),
        "core_entity_count": len({bridge["core_entity_id"] for bridge in bridges}),
        "bridges": bridges,
        "claim_boundary": "Candidate bridges are structural/source links only; they do not establish enzyme function, physiological direction, isotope tracing, or in-vivo Cannabis biosynthesis.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    return {key: result[key] for key in ("bridge_count", "target_count", "core_entity_count")}


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


def refresh_identity_set_connectivity_upstream(compounds_path: Path, identity_set_path: Path, output: Path, bq_binary: str = "bq") -> dict:
    """Extract directly characterized upstream edges for connectivity candidates.

    Candidate identity rows share only the first InChIKey block with CannabisDB
    and therefore remain a hypothesis layer. Restricting this queue to MARTS
    records keeps the first producer search tied to source-linked characterized
    enzyme records instead of importing the much larger Rhea-only neighborhood.
    """
    catalog = json.loads(compounds_path.read_text())
    compounds = catalog.get("compounds", catalog if isinstance(catalog, list) else [])
    by_prefix = {}
    exact_keys = set()
    for compound in compounds:
        key = compound.get("inchikey")
        if not key:
            continue
        by_prefix.setdefault(key.split("-")[0], []).append(compound)
        exact_keys.add(key)
    identity_report = json.loads(identity_set_path.read_text())
    prefixes = sorted({record["cannabisdb_inchikey"].split("-")[0] for record in identity_report.get("candidate_records", [])})
    if not prefixes:
        result = {"schema": "cannabis-carbon.terpene-identity-set-connectivity-upstream.v1", "edge_count": 0, "rows": [], "claim_boundary": "No candidate identity prefixes were available."}
        output.write_text(json.dumps(result, separators=(",", ":")) + "\n")
        return {"edge_count": 0, "product_count": 0}
    prefix_literals = ",".join("'" + prefix.replace("'", "\\'") + "'" for prefix in prefixes)
    exact_literals = ",".join("'" + key.replace("'", "\\'") + "'" for key in sorted(exact_keys))
    query = f"""
SELECT e.product_terpene_id, e.precursor_terpene_id, e.reaction_id,
       e.source_type, e.evidence_type, e.structure_match_mode,
       e.reaction_smarts, e.required_substrate_structures_json,
       e.missing_corpus_substrates_json, e.source_dataset, e.source_url,
       e.source_uniprot_id, e.source_genbank_id, e.source_ec_number,
       e.claim_boundary, t.inchikey AS product_inchikey, t.inchi AS product_inchi,
       t.smiles AS product_smiles, t.molecular_formula AS product_molecular_formula,
       t.carbon_count AS product_carbon_count, t.source_crossrefs AS product_source_crossrefs,
       p.inchikey AS precursor_inchikey, p.inchi AS precursor_inchi,
       p.smiles AS precursor_smiles, p.molecular_formula AS precursor_molecular_formula,
       p.carbon_count AS precursor_carbon_count, p.source_crossrefs AS precursor_source_crossrefs
FROM `terpedia-489015.terpedia_core.terpene_metabolic_map_edges_current` e
JOIN `terpedia-489015.terpedia_core.terpene_identity_set` t
  ON t.terpene_id = e.product_terpene_id
LEFT JOIN `terpedia-489015.terpedia_core.terpene_identity_set` p
  ON p.terpene_id = e.precursor_terpene_id
WHERE e.source_type = 'MARTS-DB'
  AND SUBSTR(t.inchikey, 1, 14) IN ({prefix_literals})
  AND t.inchikey NOT IN ({exact_literals})
"""
    completed = subprocess.run([bq_binary, "query", "--project_id=terpedia-489015", "--use_legacy_sql=false", "--format=prettyjson", "--max_rows=20000", query], check=True, capture_output=True, text=True)
    raw_rows = json.loads(completed.stdout)
    rows = []
    for row in raw_rows:
        prefix = (row.get("product_inchikey") or "").split("-")[0]
        candidates = by_prefix.get(prefix, [])
        rows.append({**row, "product_identity_match_status": "candidate-connectivity-inchikey", "candidate_cannabisdb_ids": sorted({c["id"] for c in candidates}), "candidate_cannabisdb_names": sorted({c.get("label") for c in candidates if c.get("label")}), "candidate_carbon_counts": sorted({c.get("carbon_atom_count") for c in candidates}), "claim_boundary": "Directly characterized Terpedia/MARTS edge joined to a connectivity-only CannabisDB identity candidate; this is a testable producer hypothesis, not exact identity or proof of endogenous Cannabis biosynthesis."})
    rows.sort(key=lambda row: (row.get("product_terpene_id") or "", row.get("reaction_id") or "", row.get("precursor_terpene_id") or ""))
    result = {"schema": "cannabis-carbon.terpene-identity-set-connectivity-upstream.v1", "source_table": EDGE_TABLE, "identity_set_table": SOURCE_TABLE, "source_identity_candidates": str(identity_set_path), "retrieved_at": date.today().isoformat(), "method": "MARTS-DB edges whose product identity shares only the first InChIKey block with a CannabisDB compound; all candidate CannabisDB IDs are retained.", "candidate_prefix_count": len(prefixes), "edge_count": len(rows), "product_count": len({row.get("product_terpene_id") for row in rows}), "reaction_count": len({row.get("reaction_id") for row in rows if row.get("reaction_id")}), "rows": rows, "claim_boundary": "These are source-linked directly characterized producer hypotheses for connectivity-only identity candidates. They do not establish stereochemical identity, reaction direction in Cannabis, enzyme specificity, endogenous biosynthesis, isotope tracing, or CO2 carbon provenance."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    return {"edge_count": result["edge_count"], "product_count": result["product_count"], "reaction_count": result["reaction_count"]}


def refresh_identity_set_candidate_expansion(connectivity_path: Path, output: Path, bq_binary: str = "bq", max_depth: int = 3) -> dict:
    """Expand connectivity-candidate precursors upstream as source-linked hypotheses."""
    seed = json.loads(connectivity_path.read_text())
    frontier = sorted({row.get("precursor_terpene_id") for row in seed.get("rows", []) if row.get("precursor_terpene_id")})
    seen_products = set(frontier)
    rows = []
    for depth in range(1, max_depth + 1):
        if not frontier:
            break
        literals = ",".join("'" + value.replace("'", "\\'") + "'" for value in frontier)
        query = f"""
SELECT e.product_terpene_id, e.precursor_terpene_id, e.reaction_id,
       e.source_type, e.evidence_type, e.structure_match_mode,
       e.reaction_smarts, e.required_substrate_structures_json,
       e.missing_corpus_substrates_json, e.source_dataset, e.source_url,
       e.source_uniprot_id, e.source_genbank_id, e.source_ec_number,
       e.claim_boundary, p.inchikey AS product_inchikey, p.inchi AS product_inchi,
       p.smiles AS product_smiles, p.molecular_formula AS product_molecular_formula,
       p.carbon_count AS product_carbon_count, p.source_crossrefs AS product_source_crossrefs,
       q.inchikey AS precursor_inchikey, q.inchi AS precursor_inchi,
       q.smiles AS precursor_smiles, q.molecular_formula AS precursor_molecular_formula,
       q.carbon_count AS precursor_carbon_count, q.source_crossrefs AS precursor_source_crossrefs
FROM `{EDGE_TABLE}` e
LEFT JOIN `{SOURCE_TABLE}` p ON p.terpene_id = e.product_terpene_id
LEFT JOIN `{SOURCE_TABLE}` q ON q.terpene_id = e.precursor_terpene_id
WHERE e.product_terpene_id IN ({literals})
"""
        completed = subprocess.run([bq_binary, "query", "--project_id=terpedia-489015", "--use_legacy_sql=false", "--format=prettyjson", "--max_rows=20000", query], check=True, capture_output=True, text=True)
        batch = json.loads(completed.stdout)
        next_frontier = set()
        for row in batch:
            product_id = row.get("product_terpene_id")
            precursor_id = row.get("precursor_terpene_id")
            if not product_id or not precursor_id:
                continue
            keep = ("product_terpene_id", "precursor_terpene_id", "reaction_id", "source_type", "evidence_type", "structure_match_mode", "reaction_smarts", "required_substrate_structures_json", "missing_corpus_substrates_json", "source_dataset", "source_url", "source_uniprot_id", "source_genbank_id", "source_ec_number", "product_inchikey", "product_inchi", "product_smiles", "product_molecular_formula", "product_carbon_count", "precursor_inchikey", "precursor_inchi", "precursor_smiles", "precursor_molecular_formula", "precursor_carbon_count")
            compact = {key: row.get(key) for key in keep if row.get(key) is not None}
            rows.append({**compact, "expansion_depth": depth, "claim_boundary": "This is a source-linked upstream expansion from a connectivity-only CannabisDB identity candidate; it is a testable hypothesis and does not establish exact stereochemical identity, Cannabis enzyme activity, endogenous biosynthesis, or CO2 provenance."})
            if precursor_id not in seen_products:
                next_frontier.add(precursor_id)
        seen_products.update(next_frontier)
        frontier = sorted(next_frontier)
    unique = {(row.get("product_terpene_id"), row.get("precursor_terpene_id"), row.get("reaction_id"), row.get("expansion_depth")): row for row in rows}
    rows = sorted(unique.values(), key=lambda row: (row.get("expansion_depth", 0), row.get("product_terpene_id") or "", row.get("reaction_id") or "", row.get("precursor_terpene_id") or ""))
    result = {"schema": "cannabis-carbon.terpene-identity-set-candidate-expansion.v1", "source_connectivity_candidates": str(connectivity_path), "source_table": EDGE_TABLE, "identity_set_table": SOURCE_TABLE, "max_depth": max_depth, "seed_product_count": len({row.get("precursor_terpene_id") for row in seed.get("rows", [])}), "expanded_edge_count": len(rows), "expanded_product_count": len({row.get("product_terpene_id") for row in rows}), "expanded_precursor_count": len({row.get("precursor_terpene_id") for row in rows}), "source_type_counts": dict(sorted(Counter(row.get("source_type") or "unknown" for row in rows).items())), "rows": rows, "claim_boundary": "Expanded edges are source-linked pathway hypotheses rooted in connectivity-only identity candidates. They remain separate from the balanced reaction inventory and directed CO2 lineage until exact identity, direction, stoichiometry, enzyme evidence, and carbon mapping are independently established."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    return {key: result[key] for key in ("seed_product_count", "expanded_edge_count", "expanded_product_count", "expanded_precursor_count")}
