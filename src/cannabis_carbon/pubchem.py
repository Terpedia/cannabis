"""Resolve CannabisDB structures against PubChem without changing identity claims."""

from __future__ import annotations

import json
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
import re
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path


PROPERTIES = "InChIKey,Title,IUPACName,CanonicalSMILES,IsomericSMILES,MolecularFormula,InChI"
BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey"
ID_EXCHANGE_URL = "https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi"


def _fetch_batch(keys: list[str], retries: int = 3) -> list[dict]:
    path = ",".join(urllib.parse.quote(key, safe="-") for key in keys)
    url = f"{BASE_URL}/{path}/property/{PROPERTIES}/JSON"
    context = ssl.create_default_context()
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, context=context, timeout=60) as response:
                payload = json.load(response)
            return payload.get("PropertyTable", {}).get("Properties", [])
        except urllib.error.HTTPError as exc:
            # PUG returns 404 when no key in a batch resolves; that is a valid
            # negative result, not permission to infer an identity.
            if exc.code == 404:
                return []
            if attempt + 1 == retries:
                raise
            time.sleep(2**attempt)
        except (urllib.error.URLError, TimeoutError):
            if attempt + 1 == retries:
                raise
            time.sleep(2**attempt)
    return []


def _fetch_bulk(keys: list[str], input_type: str = "inchikey", operator: str = "samecid", poll_seconds: int = 5, timeout_seconds: int = 1800) -> list[dict]:
    """Use PubChem's queued Identifier Exchange service for large key sets."""
    payload = urllib.parse.urlencode({
        "inputtype": input_type, "inputdsn": "", "idinput": "str",
        "idstr": "\n".join(keys), "operatortype": operator, "outputtype": "cid",
        "outputdsn": "", "method": "file-pair", "compression": "none",
        "submitjob": "Submit Job",
    }).encode()
    request = urllib.request.Request(ID_EXCHANGE_URL, data=payload, headers={"Content-Type": "application/x-www-form-urlencoded"})
    with urllib.request.urlopen(request, timeout=60) as response:
        html = response.read().decode("utf-8", "replace")
    download_match = re.search(r"https://pubchem\.ncbi\.nlm\.nih\.gov/rest/download/[^\"< ]+", html)
    reqid_match = re.search(r"reqid=(\d+)", html)
    if download_match:
        download_url = download_match.group(0).replace("&amp;", "&")
    elif reqid_match:
        download_url = None
    else:
        raise RuntimeError("PubChem Identifier Exchange did not return a request ID or download URL")
    reqid = reqid_match.group(1) if reqid_match else None
    start = re.search(r"start=([^&\"]+)", html)
    query = {"reqid": reqid, "inputtype": input_type, "inputdsn": "", "operatortype": operator, "outputtype": "cid", "outputdsn": "", "method": "file-pair", "compression": "none", "progress": "1"}
    if start:
        query["start"] = urllib.parse.unquote(start.group(1))
    if not download_url:
        status_url = f"{ID_EXCHANGE_URL}?{urllib.parse.urlencode(query)}"
        deadline = time.monotonic() + timeout_seconds
        while time.monotonic() < deadline:
            with urllib.request.urlopen(status_url, timeout=60) as response:
                status_html = response.read().decode("utf-8", "replace")
            found = re.search(r"https://pubchem\.ncbi\.nlm\.nih\.gov/rest/download/[^\"< ]+", status_html)
            if found:
                download_url = found.group(0).replace("&amp;", "&")
                break
            time.sleep(poll_seconds)
    if not download_url:
        raise TimeoutError(f"PubChem Identifier Exchange request {reqid} exceeded {timeout_seconds}s")
    with urllib.request.urlopen(download_url, timeout=60) as response:
        text = response.read().decode("utf-8", "replace")
    rows = []
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1].isdigit():
            rows.append({"query": parts[0], "CID": int(parts[1]), "InChIKey": parts[0] if input_type == "inchikey" else None})
    return rows


def resolve_pubchem(compounds_path: Path, output: Path, batch_size: int = 25, pause: float = 0.25, workers: int = 4, cache_path: Path | None = None, method: str = "batch") -> dict:
    """Resolve every CannabisDB record by exact InChIKey and retain negatives."""
    source = json.loads(compounds_path.read_text())
    compounds = source.get("compounds", source if isinstance(source, list) else [])
    records = [{"cannabisdb_id": c["id"], "inchikey": c.get("inchikey"), "smiles": c.get("smiles"), "names": sorted({x for x in [c.get("label"), *(c.get("aliases") or [])] if x}), "cannabisdb_external_ids": c.get("external_ids", {}), "cannabisdb_pubchem_cid": c.get("external_ids", {}).get("pubchem"), "status": "unresolved", "provenance": ["https://pubchem.ncbi.nlm.nih.gov/", "https://cannabisdatabase.ca/simple/download_compound_as_xml"]} for c in compounds]
    by_key = {r["inchikey"]: [] for r in records if r.get("inchikey")}
    cached = json.loads(cache_path.read_text()) if cache_path and cache_path.exists() else {}
    negative_keys = set(cached.get("negative_keys", []))
    if cached:
        for key, hits in cached.get("by_inchikey", {}).items():
            if key in by_key:
                by_key[key] = hits
    keys = sorted(key for key, hits in by_key.items() if not hits and key not in negative_keys)
    if method == "bulk":
        if keys:
            returned = _fetch_bulk(keys)
            for item in returned:
                if item.get("InChIKey") in by_key:
                    by_key[item["InChIKey"]].append(item)
            negative_keys.update(key for key in keys if not by_key[key])
        # PubChem may have standardized a structure differently (stereo,
        # charge, tautomer, or salt). Keep connectivity matches as candidates.
        fallback_records = [r for r in records if not by_key.get(r.get("inchikey"), []) and r.get("smiles")]
        by_smiles = {r["smiles"]: [] for r in fallback_records}
        if by_smiles:
            for item in _fetch_bulk(sorted(by_smiles), input_type="smiles", operator="samecon"):
                if item.get("query") in by_smiles:
                    by_smiles[item["query"]].append({"CID": item["CID"], "query": item["query"], "match_type": "same_connectivity"})
            for record in fallback_records:
                candidates = by_smiles.get(record["smiles"], [])
                if candidates:
                    record["connectivity_candidates"] = candidates
        name_records = [r for r in fallback_records if not r.get("connectivity_candidates") and r.get("names")]
        by_name = {name: [] for r in name_records for name in r["names"]}
        if by_name:
            for item in _fetch_bulk(sorted(by_name), input_type="synofiltered", operator="samecid"):
                if item.get("query") in by_name:
                    by_name[item["query"]].append({"CID": item["CID"], "query": item["query"], "match_type": "same_name"})
            for record in name_records:
                candidates = [item for name in record["names"] for item in by_name.get(name, [])]
                if candidates:
                    record["name_candidates"] = candidates
    else:
        batches = [keys[start:start + batch_size] for start in range(0, len(keys), batch_size)]
        with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
            results = executor.map(_fetch_batch, batches)
            for batch, returned in zip(batches, results):
                for item in returned:
                    key = item.get("InChIKey")
                    if key in by_key:
                        by_key[key].append(item)
                negative_keys.update(key for key in batch if not by_key[key])
                if cache_path:
                    cache_path.parent.mkdir(parents=True, exist_ok=True)
                    cache_path.write_text(json.dumps({"schema": "cannabis-carbon.pubchem-cache.v1", "by_inchikey": by_key, "negative_keys": sorted(negative_keys)}, separators=(",", ":")) + "\n")
                if pause and batch is not batches[-1]:
                    time.sleep(pause)
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps({"schema": "cannabis-carbon.pubchem-cache.v1", "by_inchikey": by_key, "negative_keys": sorted(negative_keys)}, separators=(",", ":")) + "\n")
    for record in records:
        hits = by_key.get(record.get("inchikey"), [])
        if len(hits) == 1:
            record["status"] = "resolved"
            record["pubchem"] = hits[0]
        elif len(hits) > 1:
            record["status"] = "ambiguous"
            record["pubchem_candidates"] = hits
        elif record.get("connectivity_candidates"):
            record["status"] = "candidate_connectivity"
            record["pubchem_candidates"] = record.pop("connectivity_candidates")
            record["reason"] = "same-connectivity-candidate-not-exact-identity"
        elif record.get("name_candidates"):
            record["status"] = "candidate_name"
            record["pubchem_candidates"] = record.pop("name_candidates")
            record["reason"] = "same-name-candidate-not-exact-identity"
        else:
            record["reason"] = "no-exact-inchikey-or-connectivity-match"
    summary = {
        "total": len(records),
        "resolved": sum(r["status"] == "resolved" for r in records),
        "ambiguous": sum(r["status"] == "ambiguous" for r in records),
        "candidate_connectivity": sum(r["status"] == "candidate_connectivity" for r in records),
        "candidate_name": sum(r["status"] == "candidate_name" for r in records),
        "unresolved": sum(r["status"] == "unresolved" for r in records),
        "missing_inchikey": sum(not r.get("inchikey") for r in records),
        "cannabisdb_pubchem_xref": sum(bool(r.get("cannabisdb_pubchem_cid")) for r in records),
        "xref_without_exact_resolution": sum(bool(r.get("cannabisdb_pubchem_cid")) and r["status"] not in {"resolved", "ambiguous"} for r in records),
    }
    report = {
        "schema": "cannabis-carbon.pubchem-resolution.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": str(compounds_path),
        "method": "PubChem PUG REST exact InChIKey lookup",
        "summary": summary,
        "records": records,
        "claim_boundary": "A PubChem match confirms cross-database structure identity only; it does not establish Terpedia identity, Cannabis origin, biosynthesis, enzyme function, or carbon provenance.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    return summary
