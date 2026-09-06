"""Element-level medium review; never infer nutrient need from selected fluxes."""
import hashlib
import json
from pathlib import Path


def build(inventory, evidence):
    classes={}
    for reference in evidence['references']:
        for category in ('nonmineral_essential_elements','mineral_macronutrients','mineral_micronutrients'):
            for element in reference['claims'][category]:
                if element in classes and classes[element]!=category:
                    raise ValueError('Conflicting essential element categories')
                classes[element]=category
    species=inventory['exchange_species']
    if len({c['id'] for c in species})!=len(species):
        raise ValueError('Duplicate exchange species')
    elements=set(classes)|{e for c in species for e in c['element_counts']}
    rows=[]
    for element in sorted(elements):
        containing=[c for c in species if c['element_counts'].get(element,0)>0]
        consumed=[c for c in containing if c['consuming_structure_count']>0]
        consumers=sorted({cid for c in consumed for cid in c['consuming_certificate_ids']})
        rows.append({'element':element,'reference_category':classes.get(element,'not-listed-as-essential-in-selected-reference'),
            'allowed_species_ids':[c['id'] for c in containing],
            'consumed_species_ids':[c['id'] for c in consumed],
            'consuming_certificate_ids':consumers,'consuming_structure_count':len(consumers),
            'review_status':('essential-reference-element-without-external-form' if not containing else
                'essential-reference-element-without-selected-net-demand' if not consumed else
                'essential-reference-element-with-selected-net-demand') if element in classes else
                'outside-essential-reference-review-inventory-and-exposure-scope',
            'uptake_form_supported':False,'required_amount_established':False})
    return {'schema':'cannabis-carbon.phase1-medium-elements.v1','elements':rows,
        'species':[{'compound_id':c['id'],'smiles':c['smiles'],'element_counts':c['element_counts'],
                    'role_in_saved_certificates':c['role_in_saved_certificates'],
                    'review_flags':c['review_flags']} for c in species],
        'reference_evidence':evidence,
        'summary':{'external_species':len(species),'elements_in_boundary':len({e for c in species for e in c['element_counts']}),
            'essential_reference_elements':len(classes),
            'essential_elements_without_external_form':[r['element'] for r in rows if r['review_status']=='essential-reference-element-without-external-form'],
            'essential_elements_without_selected_net_demand':[r['element'] for r in rows if r['element'] in classes and not r['consumed_species_ids']],
            'nonreference_elements_consumed':[r['element'] for r in rows if r['element'] not in classes and r['consumed_species_ids']]},
        'claim_boundary':'Descriptive element review of the selenium-forward saved certificate set. General plant-essential element evidence is not a Cannabis-specific minimum medium. Presence of an element in an allowed compound does not establish uptake, assimilation, bioavailability or the correct chemical form. Zero net consumption does not establish dispensability: catalytic pools, structural needs and biomass dilution are not modeled. Nonreference elements are not automatically toxic, nonessential in all contexts, or biosynthetic targets. No species, reactions, bounds, target denominators or pathway coverage are changed. This report does not validate the later ketone scenario or the in-progress restricted-medium scenario.'}


def run():
    paths=[Path('data/reports/phase1-medium-inventory.json'),Path('data/curation/plant-medium-evidence.json')]
    docs=[json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale evidence')
    report=build(*docs)
    report['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-medium-elements.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']))


if __name__=='__main__':
    run()
