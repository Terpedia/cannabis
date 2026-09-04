"""Refresh exact CannabisDB links from Terpedia's BigQuery identity set."""

from __future__ import annotations

import json
import subprocess
from datetime import date
from pathlib import Path


SOURCE_TABLE = "terpedia-489015.terpedia_core.terpene_identity_set"


def refresh_identity_set(compounds_path: Path, output: Path, bq_binary: str = "bq") -> dict:
    """Query the live Terpedia identity set by exact and connectivity InChIKey.

    Full-key matches are retained as exact identity evidence. Rows sharing
    only the first InChIKey block are retained in a separate candidate layer;
    that block captures connectivity but not stereochemistry, so these rows
    must never be merged into the exact layer.
    """
    catalog = json.loads(compounds_path.read_text())
    compounds = catalog.get("compounds", catalog if isinstance(catalog, list) else [])
    by_key = {}
    for compound in compounds:
        key = compound.get("inchikey")
        if key and key not in by_key:
            by_key[key] = compound
    keys = sorted(by_key)
    literals = ",".join("'" + key.replace("'", "\\'") + "'" for key in keys)
    prefixes = sorted({key.split("-")[0] for key in keys})
    query = (
        "SELECT terpene_id, identity_set_key, identity_key_type, inchi, inchikey, "
        "smiles, molecular_formula, carbon_count, source_memberships, "
        "source_record_ids, classification_statuses, classification_evidence, "
        "source_releases, manifest_uris, source_file_uris, source_crossrefs, "
        "generated_at FROM `" + SOURCE_TABLE + "` WHERE inchikey IN (" + literals + ")"
        " OR SUBSTR(inchikey, 1, 14) IN (" + ",".join("'" + prefix + "'" for prefix in prefixes) + ")"
    )
    completed = subprocess.run(
        [bq_binary, "query", "--use_legacy_sql=false", "--format=prettyjson", "--max_rows=10000", query],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = json.loads(completed.stdout)
    exact_records = []
    candidate_records = []
    exact_keys = set()
    compounds_by_prefix = {}
    for key, compound in by_key.items():
        compounds_by_prefix.setdefault(key.split("-")[0], []).append(compound)
    for row in rows:
        row_key = row.get("inchikey")
        compound = by_key.get(row_key)
        if compound is not None:
            exact_keys.add(row_key)
            exact_records.append({
                "cannabisdb_id": compound["id"],
                "cannabisdb_name": compound.get("label"),
                "cannabisdb_inchikey": compound.get("inchikey"),
                "terpedia_identity": row,
                "identity_status": "exact-inchikey",
                "match_method": "full-inchikey",
            })
            continue
        for candidate in compounds_by_prefix.get((row_key or "").split("-")[0], []):
            candidate_records.append({
                "cannabisdb_id": candidate["id"],
                "cannabisdb_name": candidate.get("label"),
                "cannabisdb_inchikey": candidate.get("inchikey"),
                "terpedia_identity": row,
                "identity_status": "candidate-connectivity-inchikey",
                "match_method": "first-inchikey-block-only",
            })
    exact_records.sort(key=lambda record: (record["cannabisdb_id"], record["terpedia_identity"].get("terpene_id", "")))
    candidate_records = [record for record in candidate_records if record["cannabisdb_inchikey"] not in exact_keys]
    candidate_records.sort(key=lambda record: (record["cannabisdb_id"], record["terpedia_identity"].get("terpene_id", "")))
    matched_candidates = {record["cannabisdb_id"] for record in candidate_records}
    report = {
        "schema": "cannabis-carbon.terpene-identity-set-match.v2",
        "source_table": SOURCE_TABLE,
        "source_project": "terpedia-489015",
        "retrieved_at": date.today().isoformat(),
        "method": "Exact full InChIKey join to CannabisDB XML-derived catalog; Terpedia identity records are identity evidence, not biosynthetic evidence.",
        "cannabisdb_compounds": len(compounds),
        "terpedia_rows_returned": len(exact_records),
        "terpedia_rows_scanned": len(rows),
        "matched_cannabisdb_compounds": len({record["cannabisdb_id"] for record in exact_records}),
        "matched_inchikeys": len({record["cannabisdb_inchikey"] for record in exact_records}),
        "connectivity_candidate_cannabisdb_compounds": len(matched_candidates),
        "connectivity_candidate_records": len(candidate_records),
        "records": exact_records,
        "candidate_records": candidate_records,
        "claim_boundary": "A verified chemical identity match does not establish endogenous Cannabis occurrence, reaction direction, enzyme function, or carbon provenance.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {key: report[key] for key in ("cannabisdb_compounds", "terpedia_rows_returned", "matched_cannabisdb_compounds", "matched_inchikeys")}
