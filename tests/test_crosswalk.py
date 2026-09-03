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
