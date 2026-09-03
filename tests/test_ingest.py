from rdkit import Chem

from cannabis_carbon.ingest import ingest_sdf


def test_ingest_marks_each_carbon_unresolved(tmp_path):
    source = tmp_path / "one.sdf"
    writer = Chem.SDWriter(str(source))
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("DATABASE_ID", "CDBTEST")
    mol.SetProp("SMILES", "CCO")
    mol.SetProp("FORMULA", "C2H6O")
    writer.write(mol)
    writer.close()
    graph = tmp_path / "graph.json"
    report = tmp_path / "report.json"
    result = ingest_sdf(source, graph, report)
    assert result["compound_count"] == 1
    assert result["carbon_atom_count"] == 2
    assert all(v == "unresolved" for v in __import__("json").loads(graph.read_text())["compounds"][0]["carbon_status"].values())
