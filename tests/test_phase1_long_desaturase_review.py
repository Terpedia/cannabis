import json
from pathlib import Path
from cannabis_carbon.phase1_long_desaturase_review import build


def test_long_desaturase_annotations_do_not_override_failed_screen(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-long-desaturase-review.json').read_text())
    assert report == build()
    assert len(report['rows']) == 2
    for row in report['rows']:
        assert row['sequence_length'] == 447
        assert row['model_eligible'] is False
        assert {'PF00173', 'PF00487'} <= {x['id'] for x in row['domain_crossrefs']}
        assert row['prior_alignments']
        assert not any(x['alignment']['passes_screen'] for x in row['prior_alignments'])
