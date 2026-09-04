"""Resolve CannabisDB structures against PubChem without changing identity claims."""

from __future__ import annotations

import json
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path


PROPERTIES = "InChIKey,Title,IUPACName,CanonicalSMILES,IsomericSMILES,MolecularFormula,InChI"
BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey"


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


def resolve_pubchem(compounds_path: Path, output: Path, batch_size: int = 25, pause: float = 0.25, workers: int = 4, cache_path: Path | None = None) -> dict:
    """Resolve every CannabisDB record by exact InChIKey and retain negatives."""
    source = json.loads(compounds_path.read_text())
    compounds = source.get("compounds", source if isinstance(source, list) else [])
    records = [{"cannabisdb_id": c["id"], "inchikey": c.get("inchikey"), "status": "unresolved", "provenance": ["https://pubchem.ncbi.nlm.nih.gov/"]} for c in compounds]
    by_key = {r["inchikey"]: [] for r in records if r.get("inchikey")}
    cached = json.loads(cache_path.read_text()) if cache_path and cache_path.exists() else {}
    negative_keys = set(cached.get("negative_keys", []))
    if cached:
        for key, hits in cached.get("by_inchikey", {}).items():
            if key in by_key:
                by_key[key] = hits
    keys = sorted(key for key, hits in by_key.items() if not hits and key not in negative_keys)
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
    for record in records:
        hits = by_key.get(record.get("inchikey"), [])
        if len(hits) == 1:
            record["status"] = "resolved"
            record["pubchem"] = hits[0]
        elif len(hits) > 1:
            record["status"] = "ambiguous"
            record["pubchem_candidates"] = hits
        else:
            record["reason"] = "no-exact-inchikey-match"
    summary = {
        "total": len(records),
        "resolved": sum(r["status"] == "resolved" for r in records),
        "ambiguous": sum(r["status"] == "ambiguous" for r in records),
        "unresolved": sum(r["status"] == "unresolved" for r in records),
        "missing_inchikey": sum(not r.get("inchikey") for r in records),
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
