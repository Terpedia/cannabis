"""Refresh exact CannabisDB links from Terpedia's BigQuery identity set."""

from __future__ import annotations

import json
import subprocess
from datetime import date
from pathlib import Path


SOURCE_TABLE = "terpedia-489015.terpedia_core.terpene_identity_set"


def refresh_identity_set(compounds_path: Path, output: Path, bq_binary: str = "bq") -> dict:
    """Query the live Terpedia identity set by exact CannabisDB InChIKey."""
    catalog = json.loads(compounds_path.read_text())
    compounds = catalog.get("compounds", catalog if isinstance(catalog, list) else [])
    by_key = {}
    for compound in compounds:
        key = compound.get("inchikey")
        if key and key not in by_key:
            by_key[key] = compound
    keys = sorted(by_key)
    literals = ",".join("'" + key.replace("'", "\\'") + "'" for key in keys)
    query = (
        "SELECT terpene_id, identity_set_key, identity_key_type, inchi, inchikey, "
        "smiles, molecular_formula, carbon_count, source_memberships, "
        "source_record_ids, classification_statuses, classification_evidence, "
        "source_releases, manifest_uris, source_file_uris, source_crossrefs, "
        "generated_at FROM `" + SOURCE_TABLE + "` WHERE inchikey IN (" + literals + ")"
    )
    completed = subprocess.run(
        [bq_binary, "query", "--use_legacy_sql=false", "--format=prettyjson", "--max_rows=10000", query],
        check=True,
        capture_output=True,
        text=True,
    )
    rows = json.loads(completed.stdout)
    records = []
    for row in rows:
        compound = by_key.get(row.get("inchikey"))
        if compound is None:
            continue
        records.append({
            "cannabisdb_id": compound["id"],
            "cannabisdb_name": compound.get("label"),
            "cannabisdb_inchikey": compound.get("inchikey"),
            "terpedia_identity": row,
        })
    records.sort(key=lambda record: record["cannabisdb_id"])
    report = {
        "schema": "cannabis-carbon.terpene-identity-set-match.v1",
        "source_table": SOURCE_TABLE,
        "source_project": "terpedia-489015",
        "retrieved_at": date.today().isoformat(),
        "method": "Exact full InChIKey join to CannabisDB XML-derived catalog; Terpedia identity records are identity evidence, not biosynthetic evidence.",
        "cannabisdb_compounds": len(compounds),
        "terpedia_rows_returned": len(records),
        "matched_cannabisdb_compounds": len(records),
        "matched_inchikeys": len({record["cannabisdb_inchikey"] for record in records}),
        "records": records,
        "claim_boundary": "A verified chemical identity match does not establish endogenous Cannabis occurrence, reaction direction, enzyme function, or carbon provenance.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {key: report[key] for key in ("cannabisdb_compounds", "terpedia_rows_returned", "matched_cannabisdb_compounds", "matched_inchikeys")}
