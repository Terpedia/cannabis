import json

from cannabis_carbon.cannabisdb_xml import enrich_compounds_with_xrefs, extract_terpedia_table


def test_xml_xrefs_are_preserved(tmp_path):
    xml = tmp_path / "compounds.xml"
    xml.write_text("<?xml version='1.0'?><compound><accession>CDB1</accession><pubchem_compound_id>123</pubchem_compound_id><chebi_id>456</chebi_id></compound><?xml version='1.0'?><compound><accession>CDB2</accession><pubchem_compound_id>Not Available</pubchem_compound_id></compound>")
    source = tmp_path / "compounds.json"
    source.write_text(json.dumps({"compounds": [{"id": "CDB1"}, {"id": "CDB2"}]}))
    out = tmp_path / "out.json"
    report = enrich_compounds_with_xrefs(xml, source, out)
    assert report["records_with_any_external_id"] == 1
    records = json.loads(out.read_text())["compounds"]
    assert records[0]["external_ids"] == {"pubchem": "123", "chebi": "456"}
    assert records[1]["external_ids"] == {}


def test_xml_table_preserves_structures_and_references(tmp_path):
    xml = tmp_path / "compounds.xml"
    xml.write_text("<compound><accession>CDB1</accession><name>Example</name><chemical_formula>C2H6O</chemical_formula><smiles>CCO</smiles><inchi>InChI=1S/C2H6O</inchi><inchikey>LFQ</inchikey><synonyms><synonym>EtOH</synonym></synonyms><general_references><reference><reference_text>Paper</reference_text><pubmed_id>123</pubmed_id></reference></general_references></compound>")
    out = tmp_path / "table.json"
    report = extract_terpedia_table(xml, out)
    assert report["compound_count"] == 1
    row = json.loads(out.read_text())["rows"][0]
    assert row["smiles"] == "CCO"
    assert row["references"][0]["pubmed_id"] == "123"
