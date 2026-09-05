"""Audit Arabidopsis purine annotations against exact Terpedia chemistry."""
import csv
import hashlib
import io
import json
import re
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from .phase1_reference_discovery import direction_families
from .phase1_scope import write_rows

QUERY = '(organism_id:3702) AND (cc_pathway:purine) AND (cc_pathway:"de novo")'
FIELDS = 'accession,reviewed,protein_name,gene_names,ec,cc_pathway,rhea'


def build(tsv, network, expanded, families):
    references = []
    for row in csv.DictReader(io.StringIO(tsv), delimiter='\t'):
        references.append({'accession':row['Entry'], 'review_status':row['Reviewed'],
            'protein_name':row['Protein names'], 'gene_names':row['Gene Names'],
            'ec_numbers':row['EC number'], 'pathway_annotation':row['Pathway'],
            'annotated_rhea_ids':re.findall(r'RHEA:\d+',row['Rhea ID']),
            'source_url':f'https://www.uniprot.org/uniprotkb/{row["Entry"]}/entry',
            'organism':'Arabidopsis thaliana', 'evidence_status':'plant-reference-annotation-not-Cannabis-activity'})
    if len({p['accession'] for p in references}) != len(references):
        raise ValueError('Duplicate reference accession')
    annotations = {m for p in references for m in p['annotated_rhea_ids']}
    audits, gaps = [], []
    for master in sorted(annotations):
        refs = [p for p in references if master in p['annotated_rhea_ids']]
        matches = []
        for r in network['reactions']:
            joins = [s for s in r['sources'] if families.get(s['source_reaction_id'],{}).get('RHEA_ID_MASTER') == master]
            if joins:
                matches.append({'reaction_id':r['id'], 'source_joins':joins,
                    'candidate_evidence_ids':expanded['candidate_reaction_evidence_ids'].get(r['id'],[])})
                if not matches[-1]['candidate_evidence_ids']:
                    participants={m['compound_id'] for side in ('left','right') for m in r[side]}
                    source_ids={s['source_reaction_id'] for s in joins}
                    gaps.append({'reaction_id':r['id'],'source_master_id':master,
                        'left':r['left'],'right':r['right'],'source_joins':joins,
                        'target_ids':[t['cannabisdb_id'] for t in network['targets'] if t['compound_id'] in participants],
                        'priority_target_ids':[],'hypothesis_ids':[],
                        'reference_matches':[{**p,'exact_source_id_matches':sorted(source_ids & set(p['annotated_rhea_ids'])),
                            'family_annotation_matches':[master], 'join_method':'explicit-published-Rhea-direction-family',
                            'direction_status':'unresolved-in-Cannabis'} for p in refs]})
        audits.append({'id':master,'reference_accessions':[p['accession'] for p in refs],
            'review_status_counts':dict(Counter(p['review_status'] for p in refs)), 'equation_matches':matches,
            'status':'missing-balanced-catalog-equation' if not matches else 'candidate-evidence-gap' if any(not m['candidate_evidence_ids'] for m in matches) else 'already-candidate-linked'})
    if len({r['reaction_id'] for r in gaps}) != len(gaps):
        raise ValueError('Ambiguous multi-family equation gap')
    return {'schema':'cannabis-carbon.phase1-plant-purine-references.v1','rows':gaps,
        'reference_proteins':references,'family_audit':audits,
        'summary':{'plant_reference_records':len(references),'review_status_counts':dict(Counter(p['review_status'] for p in references)),
            'annotated_reaction_families':len(audits),'family_status_counts':dict(Counter(a['status'] for a in audits)),
            'equation_gaps':len(gaps),'gap_reference_sequences':len({p['accession'] for r in gaps for p in r['reference_matches']}),
            'references_without_Rhea_annotation':[p['accession'] for p in references if not p['annotated_rhea_ids']]},
        'claim_boundary':'Annotation-guided plant pathway audit, not a complete pathway inventory or experimental Cannabis activity. Reviewed and unreviewed annotations remain distinct. Explicit Rhea families join exact balanced equations; no compound identities, cofactors or directions are merged. Target IDs indicate direct equation participation, not rescued routes. All input supply, compartments, direction and specificity remain unverified. Atom tracing is deferred.'}


def run():
    raw=Path('data/raw/phase1-plant-purine-references'); raw.mkdir(parents=True,exist_ok=True)
    url='https://rest.uniprot.org/uniprotkb/stream?'+urllib.parse.urlencode({'query':QUERY,'format':'tsv','fields':FIELDS})
    snapshot=raw/'arabidopsis.tsv'; metadata=raw/'request.json'
    if metadata.exists():
        lookup=json.loads(metadata.read_text())
        if lookup['url'] != url or hashlib.sha256(snapshot.read_bytes()).hexdigest()!=lookup['sha256']:
            raise ValueError('Plant reference cache mismatch')
    else:
        with urllib.request.urlopen(url,timeout=45) as response: payload=response.read()
        snapshot.write_bytes(payload)
        lookup={'url':url,'query':QUERY,'fields':FIELDS,'snapshot':str(snapshot),
            'sha256':hashlib.sha256(payload).hexdigest(),'retrieved_at':datetime.now(timezone.utc).isoformat()}
        metadata.write_text(json.dumps(lookup,indent=2)+'\n')
    paths=[Path('data/reports',n+'.json') for n in ('phase1-full-balanced-network','phase1-expanded-candidate-net')]
    directions=Path('data/raw/phase1-reference-discovery/rhea-directions.tsv')
    network,expanded=[json.loads(p.read_text()) for p in paths]
    report=build(snapshot.read_text(),network,expanded,direction_families(directions.read_text()))
    report['lookup']=lookup
    report['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths+[directions,snapshot,metadata]}
    payload=json.dumps(report,separators=(',',':'))+'\n'; sha=hashlib.sha256(payload.encode()).hexdigest()
    Path('data/reports/phase1-plant-purine-references.json').write_text(payload)
    groups=[('equation_gap','rows','reaction_id'),('reference_protein','reference_proteins','accession'),('family_audit','family_audit','id')]
    rows=[('metadata','report',{k:v for k,v in report.items() if k not in {g[1] for g in groups}})]
    for kind,key,id_key in groups: rows.extend((kind,r[id_key],r) for r in report[key])
    count=write_rows(rows,sha,Path('data/derived/phase1-plant-purine-references.ndjson'))
    print(json.dumps({'summary':report['summary'],'sha256':sha,'rows':count}))


if __name__=='__main__': run()
