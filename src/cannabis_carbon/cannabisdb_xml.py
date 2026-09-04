"""Enrich CannabisDB structure records with IDs present in its XML export."""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ID_TAGS = {
    "pubchem_compound_id": "pubchem",
    "chebi_id": "chebi",
    "kegg_id": "kegg",
    "chemspider_id": "chemspider",
    "biocyc_id": "biocyc",
    "knapsack_id": "knapsack",
    "metlin_id": "metlin",
    "foodb_id": "foodb",
    "drugbank_id": "drugbank",
    "pdb_id": "pdb",
    "bigg_id": "bigg",
    "phenol_explorer_compound_id": "phenol_explorer",
}


def _value(text: str | None) -> str | None:
    value = (text or "").strip()
    return None if not value or value.lower() in {"not available", "n/a", "na"} else value


def read_compound_xrefs(xml_path: Path) -> dict[str, dict[str, str]]:
    """Read concatenated `<compound>` XML documents keyed by CannabisDB accession."""
    text = xml_path.read_text(errors="replace")
    result = {}
    for match in re.finditer(r"<compound>.*?</compound>", text, re.DOTALL):
        element = ET.fromstring(match.group(0))
        accession = _value(element.findtext("accession"))
        if not accession:
            continue
        xrefs = {}
        for tag, prefix in ID_TAGS.items():
            value = _value(element.findtext(tag))
            if value:
                xrefs[prefix] = value
        result[accession] = xrefs
    return result


def enrich_compounds_with_xrefs(xml_path: Path, compounds_path: Path, output: Path, report_path: Path | None = None) -> dict:
    catalog = json.loads(compounds_path.read_text())
    xrefs = read_compound_xrefs(xml_path)
    counts = Counter()
    missing = []
    for compound in catalog["compounds"]:
        found = xrefs.get(compound["id"], {})
        compound["external_ids"] = {**found}
        compound["external_id_provenance"] = "CannabisDB XML export" if found else None
        for prefix in found:
            counts[prefix] += 1
        if not found:
            missing.append(compound["id"])
    result = {"schema": "cannabis-carbon.cannabisdb-xrefs.v1", "generated_at": datetime.now(timezone.utc).isoformat(), "source": str(xml_path), "compound_count": len(catalog["compounds"]), "xml_compound_count": len(xrefs), "records_with_any_external_id": len(catalog["compounds"]) - len(missing), "records_without_external_ids": len(missing), "counts_by_database": dict(sorted(counts.items())), "missing_accessions": missing, "claim_boundary": "These are identifiers supplied by CannabisDB's XML export. They support cross-database reconciliation but do not establish metabolite origin, biosynthesis, enzyme function, or carbon provenance."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(catalog, separators=(",", ":")) + "\n")
    if report_path:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(result, indent=2) + "\n")
    return result


def extract_terpedia_table(xml_path: Path, output: Path, report_path: Path | None = None) -> dict:
    """Extract a normalized, source-preserving table from the full XML export."""
    text = xml_path.read_text(errors="replace")
    rows = []
    for match in re.finditer(r"<compound>.*?</compound>", text, re.DOTALL):
        element = ET.fromstring(match.group(0))
        row = {"accession": _value(element.findtext("accession")), "name": _value(element.findtext("name")), "description": _value(element.findtext("description")), "formula": _value(element.findtext("chemical_formula")), "smiles": _value(element.findtext("smiles")), "inchi": _value(element.findtext("inchi")), "inchikey": _value(element.findtext("inchikey")), "iupac_name": _value(element.findtext("iupac_name")), "traditional_iupac": _value(element.findtext("traditional_iupac")), "molecular_weight": _value(element.findtext("average_molecular_weight")), "monoisotopic_molecular_weight": _value(element.findtext("monisotopic_molecular_weight")), "synonyms": [_value(x.text) for x in element.findall("./synonyms/synonym") if _value(x.text)], "external_ids": {}, "references": []}
        for tag, prefix in ID_TAGS.items():
            value = _value(element.findtext(tag))
            if value:
                row["external_ids"][prefix] = value
        for reference in element.findall("./general_references/reference"):
            item = {"text": _value(reference.findtext("reference_text")), "pubmed_id": _value(reference.findtext("pubmed_id"))}
            if any(item.values()):
                row["references"].append(item)
        rows.append(row)
    table = {"schema": "cannabis-carbon.terpedia-cannabisdb-table.v1", "source": str(xml_path), "compound_count": len(rows), "columns": ["accession", "name", "description", "formula", "smiles", "inchi", "inchikey", "iupac_name", "traditional_iupac", "molecular_weight", "monoisotopic_molecular_weight", "synonyms", "external_ids", "references"], "rows": rows, "claim_boundary": "This normalized table preserves CannabisDB source assertions and identifiers. It does not establish endogenous Cannabis biosynthesis, enzyme function, reaction direction, or carbon provenance."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(table, separators=(",", ":")) + "\n")
    report = {"schema": "cannabis-carbon.terpedia-cannabisdb-table-report.v1", "source": str(xml_path), "compound_count": len(rows), "records_with_external_ids": sum(bool(row["external_ids"]) for row in rows), "records_with_references": sum(bool(row["references"]) for row in rows), "claim_boundary": table["claim_boundary"]}
    if report_path:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2) + "\n")
    return report
