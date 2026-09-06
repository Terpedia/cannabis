import hashlib
import json
from pathlib import Path
from fractions import Fraction
import pytest
from cannabis_carbon.phase1_medium_inventory import build


def test_complete_saved_certificate_input_inventory():
    root=Path('data/reports'); read=lambda p:json.loads(p.read_bytes())
    report=read(root/'phase1-medium-inventory.json'); current=read(root/'phase1-selenium-forward-net.json')
    for path,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha
    layers=[read(root/n) for n in current['baseline_certificate_reports']]+[current]
    certs={c['compound_id']:c for doc in layers for c in doc.get('certificates',[])+doc.get('new_certificates',[])}
    compounds={c['id']:c for c in current['compounds']}
    rebuilt=build(current,compounds,list(certs.values()))
    assert rebuilt=={k:v for k,v in report.items() if k!='source_sha256'}
    assert report['summary']['covered_structures']==2721
    assert report['summary']['covered_target_records']==2724
    assert report['summary']['allowed_external_species']==len(current['external_exchange_compound_ids'])
    assert report['summary']['target_records_classified_as_exchange']==35
    assert report['summary']['allowed_external_species']!=35
    species={s['id']:s for s in report['exchange_species']}
    for row in report['certificate_dependencies']:
        c=certs[row['compound_id']]
        assert {p['compound_id']:Fraction(p['amount_per_target_molecule']) for p in row['net_inputs']}=={
            k:Fraction(v)/Fraction(c['target_amount']) for k,v in c['external_net_consumption'].items()}
        for p in row['net_inputs']:
            assert row['compound_id'] in species[p['compound_id']]['consuming_certificate_ids']
    assert {s['smiles'] for s in species.values() if s['carbon_count']}=={'O=C=O'}
    with pytest.raises(ValueError,match='coverage mismatch'):
        build(current,compounds,[])
    with pytest.raises(ValueError,match='Duplicate'):
        build(current,compounds,[*certs.values(),next(iter(certs.values()))])
