"""Export reviewed assignment constraints without changing reaction evidence."""
import hashlib
import json
from pathlib import Path
from .phase1_row_export import run as export_rows


def run():
    paths=[Path('data/curation/ketone-substrate-evidence-review.json'),
           Path('data/reports/phase1-ketone-stereo-hypotheses.json')]
    review,hypotheses=[json.loads(p.read_bytes()) for p in paths]
    reactions={r['id']:r for r in hypotheses['reactions']}
    compounds={c['id']:c for c in hypotheses['compounds']}
    for row in review['reviews']:
        r=reactions[row['hypothesis_id']]
        if (row['exact_target_smiles']!=compounds[row['compound_id']]['smiles'] or
            row['compound_id'] not in {p['compound_id'] for p in r['right']} or
            r['hypothesis_type']!='stereoselective-ketone-reduction' or r['enzyme_evidence_ids']):
            raise ValueError('Assignment review no longer matches exact hypothesis')
    report={**review,'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    Path('data/reports/phase1-ketone-substrate-review.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    export_rows('ketone-substrate-review')


if __name__=='__main__':
    run()
