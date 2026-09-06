import hashlib
import json
from pathlib import Path
from collections import Counter
import pytest
from rdkit import Chem
from cannabis_carbon.phase1_medium_elements import build


def test_full_element_audit_preserves_exact_species_and_unique_consumers():
    read=lambda p:json.loads(Path(p).read_bytes())
    inventory=read('data/reports/phase1-medium-inventory.json')
    evidence=read('data/curation/plant-medium-evidence.json')
    report=read('data/reports/phase1-medium-elements.json')
    assert {k:v for k,v in report.items() if k!='source_sha256'}==build(inventory,evidence)
    for p,sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
    assert len(report['species'])==102
    for c in report['species']:
        mol=Chem.AddHs(Chem.MolFromSmiles(c['smiles']))
        assert c['element_counts']==dict(Counter(a.GetSymbol() for a in mol.GetAtoms()))
    for row in report['elements']:
        inputs=[c for c in inventory['exchange_species'] if row['element'] in c['element_counts'] and c['consuming_structure_count']]
        consumers={cid for c in inputs for cid in c['consuming_certificate_ids']}
        assert set(row['consuming_certificate_ids'])==consumers
        assert row['consuming_structure_count']==len(consumers)
        assert not row['uptake_form_supported'] and not row['required_amount_established']
    assert report['summary']['essential_reference_elements']==17
    assert not report['summary']['essential_elements_without_external_form']
    assert report['summary']['essential_elements_without_selected_net_demand']==['B','Ca','K','Mg','Mn','Mo','Ni','Zn']
    assert report['summary']['nonreference_elements_consumed']==['Se']


def test_missing_element_form_is_not_replaced_with_arbitrary_species():
    evidence={'references':[{'claims':{'nonmineral_essential_elements':['C'], 'mineral_macronutrients':['K'],'mineral_micronutrients':[]}}]}
    report=build({'exchange_species':[]},evidence)
    assert report['summary']['essential_elements_without_external_form']==['C','K']
    assert not report['species']
    with pytest.raises(ValueError,match='Duplicate'):
        build({'exchange_species':[{'id':'x'},{'id':'x'}]},evidence)
