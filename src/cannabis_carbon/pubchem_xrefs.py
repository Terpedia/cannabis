"""Retrieve PubChem cross-database identifiers for CannabisDB compounds."""

from __future__ import annotations

import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


_CHEBI = re.compile(r"CHEBI:\d+", re.IGNORECASE)


def _fetch(cid: int, retries: int = 3) -> tuple[int, list[str], str | None]:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    for attempt in range(retries):
        try:
            request = Request(url, headers={"User-Agent": "Terpedia/cannabis carbon-provenance research"})
            with urlopen(request, timeout=45) as response:
                text = response.read().decode("utf-8")
            return cid, sorted({value.upper() for value in _CHEBI.findall(text)}), None
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            if attempt + 1 == retries:
                return cid, [], str(exc)
            time.sleep(1.5 * (attempt + 1))
    return cid, [], "unreachable"


def retrieve_pubchem_chebi_xrefs(pubchem_path: Path, output: Path, workers: int = 4, pause: float = 0.15) -> dict:
    """Fetch ChEBI cross-references from PubChem PUG View with provenance."""
    source = json.loads(pubchem_path.read_text())
    rows = []
    for record in source.get("records", []):
        cid = (record.get("pubchem") or {}).get("CID")
        if record.get("status") == "resolved" and cid:
            rows.append({"cannabisdb_id": record["cannabisdb_id"], "cid": int(cid), "inchikey": record.get("inchikey"), "cannabisdb_external_ids": record.get("cannabisdb_external_ids", {})})
    results = {}
    errors = {}
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = {}
        for row in rows:
            futures[pool.submit(_fetch, row["cid"])] = row
            if pause:
                time.sleep(pause / max(1, workers))
        for future in as_completed(futures):
            row = futures[future]
            cid, chebi_ids, error = future.result()
            results[cid] = chebi_ids
            if error:
                errors[str(cid)] = error
    records = []
    for row in rows:
        chebi_ids = results.get(row["cid"], [])
        records.append({**row, "pubchem_chebi_ids": chebi_ids, "new_chebi_ids": [value for value in chebi_ids if value not in {str(row["cannabisdb_external_ids"].get("chebi", "")).upper().replace("CHEBI:", "")}]})
    report = {
        "schema": "cannabis-carbon.pubchem-chebi-xrefs.v1",
        "source": str(pubchem_path),
        "method": "PubChem PUG View compound records, extracting explicit CHEBI identifiers; this is cross-reference evidence, not structure identity or biosynthesis evidence.",
        "endpoint": "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{CID}/JSON",
        "record_count": len(records),
        "records_with_chebi": sum(bool(row["pubchem_chebi_ids"]) for row in records),
        "records_with_new_chebi": sum(bool(row["new_chebi_ids"]) for row in records),
        "new_chebi_identifier_count": sum(len(row["new_chebi_ids"]) for row in records),
        "request_error_count": len(errors),
        "request_errors": errors,
        "records": records,
        "claim_boundary": "PubChem cross-references support candidate cross-database reconciliation. They do not override RDKit structure identity, establish stereochemistry, or prove endogenous Cannabis biosynthesis, pathway direction, enzyme function, or carbon provenance.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {key: report[key] for key in ("record_count", "records_with_chebi", "records_with_new_chebi", "new_chebi_identifier_count", "request_error_count")}
