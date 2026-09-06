"""Pin KEGG selenium equation and compare exact structures without identity merging."""
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from rdkit import Chem
from .phase1_balance_reference import concrete_participants
from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced

IDS = ('C00979', 'C01528', 'C05688', 'C00033')
EQUATION = 'C00979 + C01528 <=> C05688 + C00033'


def parse(reaction, molblocks):
    fields = {line[:12].strip(): line[12:].strip() for line in reaction.splitlines()
              if line[:12].strip() and len(line) > 12}
    if fields.get('EQUATION') != EQUATION:
        raise ValueError('Source equation changed; manual review required')
    structures = {}
    for key in IDS:
        mol = Chem.MolFromMolBlock(molblocks[key], sanitize=True, removeHs=True, strictParsing=True)
        if mol is None or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
            raise ValueError('Invalid concrete source structure')
        structures[key] = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    parts = concrete_participants('.'.join(structures[k] for k in IDS[:2]) + '>>' +
                                  '.'.join(structures[k] for k in IDS[2:]))
    compounds = {stable_id('structure', p['smiles']):
                 {'id': stable_id('structure', p['smiles']), **{k:v for k,v in p.items() if k != 'coefficient'}}
                 for side in parts for p in side}
    sides = [[{'compound_id': stable_id('structure', p['smiles']), 'coefficient': p['coefficient']}
              for p in side] for side in parts]
    if not balanced(sides, compounds):
        raise ValueError('Source equation is not element/isotope/charge balanced')
    return fields, structures, compounds, sides


def run():
    folder = Path('data/raw/selenium-kegg'); folder.mkdir(parents=True, exist_ok=True)
    manifest_path = folder/'retrievals.json'
    manifest = json.loads(manifest_path.read_bytes()) if manifest_path.exists() else {}
    for key in ('R03601', *IDS):
        endpoint = 'https://rest.kegg.jp/get/'+key+('/mol' if key in IDS else '')
        path = folder/(key+('.mol' if key in IDS else '.txt'))
        if key not in manifest:
            if path.exists():
                raise ValueError('Unmanifested source file')
            with urllib.request.urlopen(endpoint, timeout=45) as response:
                payload = response.read()
            path.write_bytes(payload)
            manifest[key] = {'url': endpoint, 'path': str(path),
                'retrieved_at': datetime.now(timezone.utc).isoformat(),
                'sha256': hashlib.sha256(payload).hexdigest()}
            manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
        if manifest[key]['url'] != endpoint or hashlib.sha256(path.read_bytes()).hexdigest() != manifest[key]['sha256']:
            raise ValueError('Source snapshot mismatch')
    fields, structures, compounds, sides = parse((folder/'R03601.txt').read_text(),
        {key:(folder/(key+'.mol')).read_text() for key in IDS})
    parent_path = Path('data/reports/phase1-local-speciation-net.json')
    parent = json.loads(parent_path.read_bytes())
    paths = [Path('data/reports/phase1-full-balanced-network.json'), Path('data/reports/phase1-marts-completions.json'), parent_path] + [Path('data/reports')/p for p in parent['baseline_certificate_reports']]
    current = {}
    for path in paths:
        doc = json.loads(path.read_bytes())
        for source, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(source).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale model source')
        for c in doc.get('compounds', []):
            if c['id'] in current and current[c['id']]['smiles'] != c['smiles']:
                raise ValueError('Exact structure conflict')
            current[c['id']] = c
    endpoints = [{'kegg_id':key, 'source_url':manifest[key]['url'], 'smiles':smiles,
        'compound_id':stable_id('structure', smiles),
        'exact_model_match':stable_id('structure', smiles) in current,
        'target_records':[t for t in parent['targets'] if t['compound_id']==stable_id('structure', smiles)]}
        for key,smiles in structures.items()]
    report = {'schema':'cannabis-carbon.phase1-selenium-source-identity.v1',
        'source_url':'https://www.kegg.jp/entry/R03601', 'source_fields':fields,
        'source_name_definition_conflict':'NAME mentions hydrogen sulfide; DEFINITION and exact EQUATION specify hydrogen selenide. Retain conflict pending curation.',
        'endpoints':endpoints, 'compounds':list(compounds.values()),
        'reaction':{'left':sides[0], 'right':sides[1], 'source_ec':fields.get('ENZYME'),
                    'balance_status':'independently-element-isotope-charge-balanced',
                    'direction_status':'source-reversible-not-Cannabis-physiological-direction'},
        'retrievals':manifest,
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary':{'source_endpoints':len(endpoints), 'exact_model_matches':sum(e['exact_model_match'] for e in endpoints),
                   'coverage_gain_claimed':0},
        'claim_boundary':'Source identity and balance audit only. No names-only joins, protonation merges, model insertion, enzyme assignment or CO2 pathway gain. KEGG records a general reaction, not Cannabis activity. Charged model forms require explicit separately balanced speciation hypotheses.'}
    Path('data/reports/phase1-selenium-source-identity.json').write_text(json.dumps(report, separators=(',',':'))+'\n')
    print(json.dumps(report['summary']), flush=True)
    print(json.dumps(endpoints), flush=True)


if __name__ == '__main__':
    run()
