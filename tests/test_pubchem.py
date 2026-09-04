import json

from cannabis_carbon import pubchem


def test_pubchem_resolver_retains_negative_records(tmp_path, monkeypatch):
    source = tmp_path / "compounds.json"
    source.write_text(json.dumps({"compounds": [
        {"id": "CDB1", "inchikey": "AAA-BBB-CC"},
        {"id": "CDB2", "inchikey": "DDD-EEE-FF"},
        {"id": "CDB3"},
    ]}))
    monkeypatch.setattr(pubchem, "_fetch_batch", lambda keys: [{"CID": 1, "InChIKey": "AAA-BBB-CC", "Title": "x"}])
    out = tmp_path / "pubchem.json"
    summary = pubchem.resolve_pubchem(source, out, batch_size=100, pause=0)
    assert summary == {"total": 3, "resolved": 1, "ambiguous": 0, "candidate_connectivity": 0, "unresolved": 2, "missing_inchikey": 1}
    records = json.loads(out.read_text())["records"]
    assert records[0]["pubchem"]["CID"] == 1
    assert records[1]["reason"] == "no-exact-inchikey-or-connectivity-match"
