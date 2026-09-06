"""Preserve CannabisDB name/structure assertions for unresolved PG targets."""
import gzip
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from .phase1_lipid_acylation import canonical
from .phase1_acyl_input_identity_review import alkene_geometry

BOUNDARY = ('Source-assertion audit only. Names/descriptions and encoded structures '
    'are separate evidence; naming E/Z geometry does not silently repair SMILES, '
    'InChI or historical target identity. Expected occurrence is not observed '
    'Cannabis metabolomics evidence. No reactions, identities, source records or '
    'coverage counts are changed. Any future curated identity must retain both '
    'original assertions and explicit resolution provenance.')


def build(review, current, xml_text):
    requested = {cid for row in review['inputs'] for cid in row['affected_cannabisdb_ids']}
    targets = {t['cannabisdb_id']:t for t in current['targets']}
    cs = {c['id']:c for c in current['compounds']}
    rows = []
    for match in re.finditer(r'<compound>.*?</compound>',xml_text,re.S):
        e = ET.fromstring(match.group())
        accession = e.findtext('accession')
        if accession not in requested:
            continue
        source = {k:(e.findtext(k) or '').strip() for k in ('accession','name','smiles','inchi','inchikey',
            'iupac_name','traditional_iupac','description','pubchem_compound_id','chebi_id')}
        source['references'] = [{k:(ref.findtext(k) or '').strip() for k in ('reference_text','pubmed_id')}
                                for ref in e.findall('./general_references/reference')]
        mol = Chem.MolFromSmiles(source['smiles'])
        inchi_mol = Chem.MolFromInchi(source['inchi']) if source['inchi'] else None
        if mol is None:
            raise ValueError('Unparseable source SMILES')
        encoded = alkene_geometry(mol)
        geometry_tokens = re.findall(r'\d+[EZ]',source['name'])
        absent = sum(g['stereo']=='STEREONONE' for g in encoded)
        target = targets[accession]
        rows.append({'cannabisdb_id':accession,'compound_id':target['compound_id'],
            'source_url':'https://cannabisdatabase.ca/compounds/'+accession,
            'source_fields':source,'xml_fragment_sha256':hashlib.sha256(match.group().encode()).hexdigest(),
            'xml_smiles_matches_retained_structure':canonical(mol)==cs[target['compound_id']]['smiles'],
            'source_name_geometry_tokens':geometry_tokens,'smiles_alkene_geometry':encoded,
            'unassigned_smiles_alkene_bonds':absent,
            'inchi_matches_source_smiles':canonical(inchi_mol)==canonical(mol) if inchi_mol else None,
            'name_structure_status':'name-specifies-geometry-omitted-from-structure' if geometry_tokens and absent else 'requires-manual-review',
            'occurrence_assertion':'source-says-expected-in-Cannabis-not-detection' if 'expected to be in cannabis' in source['description'].lower() else 'not-classified-by-this-audit',
            'claim_boundary':BOUNDARY})
    if {r['cannabisdb_id'] for r in rows} != requested or len(rows)!=len(requested):
        raise ValueError('Incomplete or duplicate accession extraction')
    rows.sort(key=lambda r:r['cannabisdb_id'])
    return {'schema':'cannabis-carbon.phase1-pg-source-identity-audit.v1','targets':rows,'claim_boundary':BOUNDARY,
        'summary':{'requested_records':len(requested),'source_records':len(rows),
            'name_structure_status_counts':dict(Counter(r['name_structure_status'] for r in rows)),
            'occurrence_assertion_counts':dict(Counter(r['occurrence_assertion'] for r in rows)),
            'xml_smiles_matching_retained_structure':sum(r['xml_smiles_matches_retained_structure'] for r in rows),
            'inchi_matching_smiles':sum(r['inchi_matches_source_smiles'] is True for r in rows),
            'new_reactions':0,'new_CO2_route_claims':0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-acyl-input-identity-review.json'),
             Path('data/reports/phase1-glycerophospholipid-net.json'),
             Path('data/terpedia/cannabisdb-compounds.xml.gz')]
    report = build(json.loads(paths[0].read_bytes()),json.loads(paths[1].read_bytes()),
                   gzip.decompress(paths[2].read_bytes()).decode())
    report['source_sha256'] = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-pg-source-identity-audit.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)


if __name__ == '__main__':
    run()
