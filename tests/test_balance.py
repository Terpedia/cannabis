import gzip
import json

from cannabis_carbon.balance import audit_balances


def test_phase1_requires_both_element_and_charge_balance(tmp_path):
    network = {"entities": [{"id": "r:1", "type": "biochemical_reaction", "attributes": {"elementBalance": {"status": "balanced"}, "chargeBalance": {"status": "not_auditable"}}}], "statements": []}
    source = tmp_path / "n.json.gz"
    with gzip.open(source, "wt") as h: json.dump(network, h)
    out = tmp_path / "report.json"
    summary = audit_balances(source, out)
    assert summary == {"fully_balanced": 0, "imbalanced": 0, "not_auditable": 1, "element_balanced": 1, "charge_balanced": 0}
