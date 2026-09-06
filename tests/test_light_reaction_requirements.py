import json
from pathlib import Path
from rdkit import Chem


def test_exact_light_source_joins_and_reactive_part_context():
    root=Path('.')
    annotations=json.loads((root/'data/curation/light-reaction-requirements.json').read_bytes())
    network=json.loads((root/'data/reports/phase1-full-balanced-network.json').read_bytes())
    reactions={r['id']:r for r in network['reactions']}; compounds={c['id']:c for c in network['compounds']}
    assert len(annotations['reactions'])==2
    for a in annotations['reactions']:
        r=reactions[a['model_reaction_id']]
        matches=[s for s in r['sources'] if s['source_reaction_id'].upper()==a['source_forward_id']]
        assert matches
        assert {s['source_left_corresponds_to'] for s in matches}=={'right'}
        assert a['source_forward_model_direction']=='hypothetical-right-to-left'
        assert len(a['energy_inputs'])==1
        energy=a['energy_inputs'][0]
        assert energy['entity_type']=='photon' and energy['source_id']=='CHEBI:30212'
        assert energy['coefficient']==(1 if a['source_master_id']=='RHEA:30407' else 4)
        species={compounds[p['compound_id']]['smiles'] for side in ('left','right') for p in r[side]}
        for c in a['reactive_part_context']:
            assert Chem.MolToSmiles(Chem.MolFromSmiles(c['reactive_part_smiles'])) in species
        assert '*' not in species
