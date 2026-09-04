import gzip
import json
from rdkit import Chem

from cannabis_carbon.crosswalk import build_crosswalk


def test_crosswalk_matches_exact_inchikey(tmp_path):
    sdf = tmp_path / "cdb.sdf"
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("DATABASE_ID", "CDB1")
    mol.SetProp("FORMULA", "C2H6O")
    Chem.SDWriter(str(sdf)).write(mol)
    network = {"entities": [{"id": "chebi:1", "type": "metabolite", "label": "ethanol", "attributes": {"canonicalSmiles": "CCO"}}], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    output = tmp_path / "crosswalk.json"
    result = build_crosswalk(sdf, network_path, output)
    assert result["exact_matches"] == 1
    assert result["unmatched"] == 0
    assert result["cannabisdb_unmatched"] == 0
    assert result["terpedia_unmatched"] == 0


def test_crosswalk_retains_connectivity_candidate_separately(tmp_path):
    sdf = tmp_path / "cdb.sdf"
    mol = Chem.MolFromSmiles("C[C@H](O)F")
    mol.SetProp("DATABASE_ID", "CDB1")
    Chem.SDWriter(str(sdf)).write(mol)
    network = {"entities": [{"id": "chebi:1", "type": "metabolite", "label": "fluoroethanol", "attributes": {"canonicalSmiles": "CC(O)F"}}], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    result = build_crosswalk(sdf, network_path, tmp_path / "out.json")
    assert result["exact_matches"] == 0
    assert result["connectivity_candidate_matches"] == 1


def test_crosswalk_retains_tautomer_candidate_separately(tmp_path):
    sdf = tmp_path / "cdb.sdf"
    mol = Chem.MolFromSmiles("CC(O)=C")
    mol.SetProp("DATABASE_ID", "CDB1")
    Chem.SDWriter(str(sdf)).write(mol)
    network = {"entities": [{"id": "chebi:1", "type": "metabolite", "label": "acetone tautomer", "attributes": {"canonicalSmiles": "CC(=O)C"}}], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    result = build_crosswalk(sdf, network_path, tmp_path / "out.json")
    assert result["exact_matches"] == 0
    assert result["tautomer_candidate_matches"] == 1


def test_crosswalk_retains_formula_compatible_name_candidate_separately(tmp_path):
    sdf = tmp_path / "cdb.sdf"
    mol = Chem.MolFromSmiles("CCC=O")
    mol.SetProp("DATABASE_ID", "CDB1")
    mol.SetProp("FORMULA", "C3H6O")
    mol.SetProp("GENERIC_NAME", "acetone")
    Chem.SDWriter(str(sdf)).write(mol)
    network = {"entities": [{"id": "chebi:1", "type": "metabolite", "label": "acetone", "attributes": {"canonicalSmiles": "CC(=O)C", "molecularFormula": "C3H6O"}}], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    result = build_crosswalk(sdf, network_path, tmp_path / "out.json")
    assert result["exact_matches"] == 0
    assert result["name_candidate_matches"] == 1


def test_crosswalk_includes_source_linked_reaction_addition_entities(tmp_path):
    sdf = tmp_path / "cdb.sdf"
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("DATABASE_ID", "CDB1")
    Chem.SDWriter(str(sdf)).write(mol)
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle:
        json.dump({"entities": [], "statements": []}, handle)
    (tmp_path / "reaction-additions.json").write_text(json.dumps({"entities": [
        {"id": "cannabisdb:CDB1", "type": "metabolite", "label": "ethanol", "attributes": {"canonicalSmiles": "CCO"}}
    ], "statements": []}))
    result = build_crosswalk(sdf, network_path, tmp_path / "out.json")
    assert result["exact_matches"] == 1
