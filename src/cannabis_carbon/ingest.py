from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem


def ingest_sdf(source: Path, graph_out: Path, report_out: Path) -> dict:
    """Convert a CannabisDB SDF into source-linked compound nodes and a gap report."""
    compounds = []
    invalid = []
    for index, mol in enumerate(Chem.SDMolSupplier(str(source), removeHs=False, sanitize=True)):
        if mol is None:
            invalid.append(index)
            continue
        props = {name: mol.GetProp(name) for name in mol.GetPropNames()}
        compound_id = props.get("DATABASE_ID", f"sdf-record:{index + 1}")
        carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
        compounds.append({
            "id": compound_id,
            "label": compound_id,
            "smiles": props.get("SMILES", Chem.MolToSmiles(mol)),
            "inchikey": props.get("INCHI_KEY"),
            "formula": props.get("FORMULA"),
            "carbon_atom_count": len(carbon_atoms),
            "carbon_status": {str(atom): "unresolved" for atom in carbon_atoms},
            "source": "Cannabis Compound Database",
            "source_url": "https://cannabisdatabase.ca",
        })
    graph_out.parent.mkdir(parents=True, exist_ok=True)
    report_out.parent.mkdir(parents=True, exist_ok=True)
    graph_out.write_text(json.dumps({"schema": "cannabis-carbon.compounds.v1", "compounds": compounds}, separators=(",", ":")) + "\n")
    report = {
        "schema": "cannabis-carbon.coverage.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_file": str(source),
        "compound_count": len(compounds),
        "invalid_records": invalid,
        "carbon_atom_count": sum(item["carbon_atom_count"] for item in compounds),
        "status_counts": {"supported": 0, "candidate": 0, "inferred": 0, "unresolved": sum(item["carbon_atom_count"] for item in compounds)},
        "claim_boundary": "A structure catalog does not establish biosynthesis or carbon origin.",
    }
    report_out.write_text(json.dumps(report, indent=2) + "\n")
    return report
