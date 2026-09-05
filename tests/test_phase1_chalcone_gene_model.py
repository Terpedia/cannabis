import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_chalcone_gene_model import parse_record


def test_chalcone_source_prediction_preserves_exact_sequence_and_coding_order(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-chalcone-gene-model.json').read_text())
    payload = Path('data/raw/chalcone-gene-model/KAF4401769.1.gb').read_bytes()
    parsed = parse_record(payload)
    assert report['source_record'] == parsed
    annotation = json.loads(Path('data/raw/chalcone-annotations/A0A7J6I409.json').read_text())
    assert parsed['sequence'] == annotation['sequence']['value']
    assert len(parsed['sequence']) == 346
    assert parsed['coding_length_nt'] == 1041
    exons = parsed['genomic_order_exons']
    assert len(exons) == 6
    assert {(e['accession'], e['strand']) for e in exons} == {('JAATIQ010000011.1', -1)}
    assert [(e['start'], e['end']) for e in exons] == [(7834077,7834385),(7834492,7834715),(7834799,7834957),(7835397,7835470),(7835726,7835844),(7836053,7836208)]
    assert parsed['transcript_order_exons'] == list(reversed(exons))
    assert report['model_eligible'] is False
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    with pytest.raises(ValueError):
        parse_record(payload.replace(b'complement(join(', b'join('))
