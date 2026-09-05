"""Versioned union of catalog and reference-backfill homology evidence."""
import hashlib
import json
from pathlib import Path
from .phase1_catalog_evidence import assemble
from .phase1_screened_overlay import build_overlay
from .phase1_scope import write_rows


def build(catalog, previous, search):
    added = build_overlay({'reactions':catalog['reactions'], 'hypotheses':[]}, search,
        'https://github.com/Terpedia/cannabis/blob/main/data/reports/phase1-backfill-protein-search.json')
    for e in added['enzyme_evidence']:
        e['evidence_class'] = 'reference-backfill-direction-unresolved-homology-candidate'
    combined = previous['enzyme_evidence'] + added['enzyme_evidence']
    output = assemble(catalog, {'enzyme_evidence':combined, 'summary':{
        'distinct_cannabis_proteins':len({p['accession'] for e in combined for p in e['screened_proteins']})}})
    output['supplement_version'] = 'catalog-and-reference-backfill'
    output['summary']['backfill_candidate_equations'] = len(added['enzyme_evidence'])
    output['summary']['previous_candidate_equations'] = len(previous['enzyme_evidence'])
    output['comparison_boundary'] = 'Union relative to the unchanged original catalog diagnostic, not a second application on top of the first supplement. All previous evidence records are preserved verbatim. The 19 backfill additions do not fully candidate-link another selected net certificate.'
    return output


def run():
    paths = [Path('data/reports', n+'.json') for n in ('phase1-catalog-net-gaps',
        'phase1-catalog-evidence', 'phase1-backfill-protein-search', 'phase1-reference-backfill')]
    catalog, previous, search, references = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for report in (previous,references):
        for path,digest in report['source_sha256'].items():
            if path in hashes and hashes[path] != digest:
                raise ValueError('Combined evidence source mismatch')
    if search['source_discovery_sha256'] != hashes[str(paths[3])]:
        raise ValueError('Backfill search discovery mismatch')
    output = build(catalog,previous,search); output['source_sha256'] = hashes
    payload = json.dumps(output,separators=(',',':'))+'\n'
    Path('data/reports/phase1-combined-catalog-evidence.json').write_text(payload)
    groups = [('enzyme_evidence','enzyme_evidence','id'), ('certificate_update','certificate_updates','compound_id'),
        ('target_update','target_updates','cannabisdb_id')]
    rows = [('metadata','supplement',{k:v for k,v in output.items() if k not in {g[1] for g in groups}})]
    for kind,collection,key in groups:
        rows.extend((kind,r[key],r) for r in output[collection])
    sha = hashlib.sha256(payload.encode()).hexdigest()
    count = write_rows(rows,sha,Path('data/derived/phase1-combined-catalog-evidence.ndjson'))
    print(json.dumps({'summary':output['summary'],'sha256':sha,'bytes':len(payload.encode()),'rows':count}))


if __name__ == '__main__':
    run()
