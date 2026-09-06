"""Balanced chemistry for separately named PG alternatives, never target replacements."""
import hashlib
import json
from pathlib import Path
from rdkit import Chem, RDLogger
from .phase1_glycerophospholipid_speciation import build as speciation
from .phase1_glycerophospholipid_synthesis import build as synthesis
from .phase1_pg_named_alternatives import BOUNDARY


def build(current, alternatives, catalog):
    compounds = {c['id']:dict(c) for c in current['compounds']}
    probes = []
    for row in alternatives['targets']:
        cid = row['alternative_compound_id']; smiles = row['alternative_smiles']
        if cid in compounds and compounds[cid]['smiles']!=smiles:
            raise ValueError('Identity conflict')
        compounds.setdefault(cid,{'id':cid,'smiles':smiles,'formula':row['formula'],
            'formal_charge':row['formal_charge'],'carbon_count':sum(a.GetAtomicNum()==6 for a in Chem.MolFromSmiles(smiles).GetAtoms())})
        probes.append({'cannabisdb_id':row['cannabisdb_id'],'label':row['source_name'],
            'compound_id':cid,'net_status':'name-derived-alternative-not-yet-tested'})
    probe_inventory = {**current,'compounds':list(compounds.values()),'targets':probes}
    audit = speciation(probe_inventory,catalog)
    if {t['compound_id'] for t in audit['targets']}!={t['compound_id'] for t in probes}:
        raise ValueError('Not every named alternative has an exact source scaffold')
    report = synthesis(current,audit,catalog)
    report['schema'] = 'cannabis-carbon.phase1-pg-named-reactions.v1'
    report['claim_boundary'] = BOUNDARY + ' ' + report['claim_boundary']
    report['identity_alternatives'] = alternatives['targets']
    report['speciation_audit'] = audit
    return report


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-glycerophospholipid-net.json'),
             Path('data/reports/phase1-pg-named-alternatives.json'),
             Path('data/raw/phase1-balance-reference-catalog.json')]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs[:2]:
        for p,sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Changed source snapshot')
    report = build(*docs)
    report['source_sha256'] = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-pg-named-reactions.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)


if __name__ == '__main__':
    run()
