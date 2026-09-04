import gzip
import json

from cannabis_carbon.balance import audit_balances


def test_phase1_requires_both_element_and_charge_balance(tmp_path):
    network = {"entities": [{"id": "r:1", "type": "biochemical_reaction", "attributes": {"elementBalance": {"status": "balanced"}, "chargeBalance": {"status": "not_auditable"}}}], "statements": []}
    source = tmp_path / "n.json.gz"
    with gzip.open(source, "wt") as h: json.dump(network, h)
    out = tmp_path / "report.json"
    summary = audit_balances(source, out)
    assert summary == {"fully_balanced": 0, "imbalanced": 0, "not_auditable": 1, "element_balanced": 1, "charge_balanced": 0, "computed_fully_balanced": 0, "computed_imbalanced": 0, "reaction_smiles_fully_balanced": 0, "reaction_smiles_imbalanced": 0}


def test_generic_rhea_substituents_are_not_auditable(tmp_path):
    network = {
        "entities": [
            {"id": "r:1", "type": "biochemical_reaction", "attributes": {}},
            {"id": "c:1", "type": "metabolite", "attributes": {"molecularFormula": "C5H5O8PR3", "formalCharge": -1, "canonicalSmiles": "[1*]C(=O)O"}},
            {"id": "c:2", "type": "metabolite", "attributes": {"molecularFormula": "C1H2O1", "formalCharge": 0, "canonicalSmiles": "CO"}},
        ],
        "statements": [
            {"subjectId": "r:1", "predicate": "has_reactant", "objectEntityId": "c:1"},
            {"subjectId": "r:1", "predicate": "has_product", "objectEntityId": "c:2"},
        ],
    }
    source = tmp_path / "n.json.gz"
    with gzip.open(source, "wt") as handle:
        json.dump(network, handle)
    summary = audit_balances(source, tmp_path / "report.json")
    assert summary["imbalanced"] == 0
    assert summary["not_auditable"] == 1


def test_concrete_reaction_smiles_can_close_formula_gap(tmp_path):
    network = {
        "entities": [
            {"id": "r:1", "type": "biochemical_reaction", "attributes": {"reactionSmiles": "C=C.O>>CCO"}},
            {"id": "c:1", "type": "metabolite", "attributes": {}},
            {"id": "c:2", "type": "metabolite", "attributes": {}},
            {"id": "c:3", "type": "metabolite", "attributes": {}},
        ],
        "statements": [
            {"subjectId": "r:1", "predicate": "has_reactant", "objectEntityId": "c:1"},
            {"subjectId": "r:1", "predicate": "has_reactant", "objectEntityId": "c:2"},
            {"subjectId": "r:1", "predicate": "has_product", "objectEntityId": "c:3"},
        ],
    }
    source = tmp_path / "n.json.gz"
    with gzip.open(source, "wt") as handle:
        json.dump(network, handle)
    summary = audit_balances(source, tmp_path / "report.json")
    assert summary["fully_balanced"] == 1
    assert summary["reaction_smiles_fully_balanced"] == 1
