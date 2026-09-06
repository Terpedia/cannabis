"""Exact exchange dependencies of saved certificates, not a minimum-medium claim."""
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from rdkit import Chem, RDLogger
from .phase1_glycerophospholipid_input_audit import assemble_current
from .phase1_reaction_completion_net import validate_certificate


def build(current, compounds, certificates):
    exchange=set(current['external_exchange_compound_ids'])
    by_id={c['compound_id']:c for c in certificates}
    if len(by_id)!=len(certificates):
        raise ValueError('Duplicate certificate identity')
    targets=defaultdict(list)
    for t in current['targets']:
        if t['net_status']=='exact-net-conversion-hypothesis':
            targets[t['compound_id']].append(t['cannabisdb_id'])
    if set(targets)!=set(by_id):
        raise ValueError('Certificate coverage mismatch')
    consumed=defaultdict(list); exported=defaultdict(list)
    rows=[]
    for cid, cert in sorted(by_id.items()):
        inputs=[]; outputs=[]
        for field,index,destination in [('external_net_consumption',consumed,inputs),('net_exports',exported,outputs)]:
            for species,amount in cert[field].items():
                if species not in exchange:
                    if field=='external_net_consumption':
                        raise ValueError('Nonexchange input')
                    continue
                if Fraction(amount)<=0:
                    raise ValueError('Nonpositive exchange amount')
                normalized=str(Fraction(amount)/Fraction(cert['target_amount']))
                destination.append({'compound_id':species,'amount_per_target_molecule':normalized})
                index[species].append(cid)
        rows.append({'compound_id':cid,'cannabisdb_ids':targets[cid],
                     'net_inputs':inputs,'net_external_outputs':outputs})
    species=[]
    for cid in sorted(exchange):
        c=compounds[cid]; mol=Chem.MolFromSmiles(c['smiles'])
        elements=Counter(a.GetSymbol() for a in Chem.AddHs(mol).GetAtoms())
        flags=[]
        if elements.get('Fe',0)>1 and elements.get('S',0)>1:
            flags.append('iron-sulfur-cluster-exchange-review-internal-carrier-synthesis')
        if c['smiles'] in ('OO','[O][O-]'):
            flags.append('reactive-oxygen-exchange-review-energy-redox-boundary')
        if c['smiles']=='[O-]P([O-])([O-])=[Se]':
            flags.append('activated-selenium-donor-review-internal-synthesis')
        if c['smiles']=='[H+]':
            flags.append('explicit-proton-balance-not-a-specified-buffer-recipe')
        species.append({**c,'element_counts':dict(elements),'review_flags':flags,
            'consuming_certificate_ids':consumed[cid],
            'consuming_structure_count':len(consumed[cid]),
            'consuming_record_count':sum(len(targets[t]) for t in consumed[cid]),
            'exporting_structure_count':len(exported[cid]),
            'role_in_saved_certificates':'net-input' if consumed[cid] else 'export-only' if exported[cid] else 'unused',
            'nutrient_essentiality':'not-established-by-this-audit'})
    return {'schema':'cannabis-carbon.phase1-medium-inventory.v1','exchange_species':species,
        'certificate_dependencies':rows,
        'summary':{'inventory_target_records':len(current['targets']),
            'covered_target_records':sum(map(len,targets.values())), 'covered_structures':len(by_id),
            'allowed_external_species':len(exchange),
            'allowed_noncarbon_external_species':sum(not c['carbon_count'] for c in species),
            'target_records_classified_as_exchange':sum(t['net_status'].startswith('explicit-exchange-species') for t in current['targets']),
            'used_net_input_species':sum(bool(consumed[c]) for c in exchange),
            'used_noncarbon_net_input_species':sum(bool(consumed[c['id']]) and not c['carbon_count'] for c in species),
            'flagged_used_input_species':sum(bool(c['review_flags']) and c['consuming_structure_count']>0 for c in species)},
        'claim_boundary':'Union of net inputs used by one saved certificate per covered exact structure. This is an upper-bound input list for those specific certificates, not a minimum medium, essential nutrient list or growth recipe. Alternative pathways may use fewer or different inputs. Certificates permit regenerated pre-existing pools; catalytic minerals with zero net consumption, biomass growth, transport, light, energy, concentrations, pH and counterions are not established. The full inventory includes environmental/exposure compounds and cannot automatically define a normal plant growth medium. No exchange bounds or pathway coverage were changed.',
        'next_steps':['Curate admissible external nutrient forms separately from internal carriers and reactive intermediates.',
            'Split synthesis targets from uptake/exposure records while preserving the full inventory denominator.',
            'Optimize jointly over allowed input species and full reaction fluxes, preserving all required substrates and target production constraints.',
            'Prove minimality or label the result inclusion-minimal/upper-bound; replay every target under the selected medium.',
            'Model catalytic mineral demand, biomass, light/energy, transport and starting pools before calling the result a plant minimum defined medium.']}


def run():
    RDLogger.DisableLog('rdApp.warning')
    root=Path('data/reports'); read=lambda p:json.loads(p.read_bytes())
    current_path=root/'phase1-selenium-forward-net.json'; current=read(current_path)
    paths=[root/'phase1-full-balanced-network.json',root/'phase1-marts-completions.json',current_path]+[root/n for n in current['baseline_certificate_reports']]
    docs=[read(p) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale source')
    _,compounds,model=assemble_current(docs[0],docs[1],current,docs[3:])
    certs={c['compound_id']:c for doc in [*docs[3:],current] for c in doc.get('certificates',[])+doc.get('new_certificates',[])}
    steps={s['id']:s for s in model.steps}
    for c in certs.values():
        validate_certificate(c,steps,compounds,set(current['external_exchange_compound_ids']),current['co2_compound_id'])
    report=build(current,compounds,list(certs.values()))
    report['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (root/'phase1-medium-inventory.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)
    print(json.dumps([{k:c[k] for k in ('smiles','consuming_structure_count','review_flags')} for c in report['exchange_species'] if c['consuming_structure_count']]),flush=True)


if __name__=='__main__':
    run()
