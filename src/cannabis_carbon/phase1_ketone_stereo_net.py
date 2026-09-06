"""Full-inventory forward-only sensitivity for explicit ketone redox hypotheses."""
import hashlib
import json
from pathlib import Path
from .phase1_lipid_acylation_net import build


def run():
    root=Path('data/reports'); read=lambda p:json.loads(p.read_bytes())
    parent_path=root/'phase1-selenium-forward-net.json'; parent=read(parent_path)
    paths=[root/'phase1-full-balanced-network.json',root/'phase1-marts-completions.json',parent_path,
           root/'phase1-ketone-stereo-hypotheses.json']+[root/n for n in parent['baseline_certificate_reports']]
    docs=[read(p) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=sha:
                raise ValueError('Stale model source')
    previous=docs[4:]
    certs={c['compound_id']:c for doc in previous for c in doc.get('certificates',[])+doc.get('new_certificates',[])}
    if len(certs)+len(parent['new_certificates'])!=2721:
        raise ValueError('Missing historical witnesses')
    original={'certificates':list(certs.values()),'external_exchange_compound_ids':parent['external_exchange_compound_ids']}
    layers=[{**doc,'forbidden_step_ids':doc.get('forbidden_step_ids',[])} for doc in [*previous,parent] if 'added_reactions' in doc]
    report=build(docs[0],docs[1],original,parent,docs[3],prior_layers=layers,
                 extra_forward_types=('secondary-alcohol-oxidation','stereoselective-ketone-reduction'))
    report.update({'schema':'cannabis-carbon.phase1-ketone-stereo-net.v1',
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'baseline_certificate_reports':parent['baseline_certificate_reports']+[parent_path.name],
        'lipid_evidence_report':'phase1-ketone-stereo-hypotheses.json',
        'claim_boundary':'Reaction-class-analogy sensitivity, not substrate-specific experimental evidence or demonstrated Cannabis activity. Each newly added oxidation/reduction equation runs only in its proposed forward direction; exact existing-equation joins keep inherited directions. No direct stereoisomerization edge, identity merge or external input change. Positive results are conditional net-conversion certificates under the inherited permissive inorganic boundary, allowing regenerated pre-existing pools. They do not establish a minimum medium, startup, energetics, light, compartments, transport or enzyme stereoselectivity.'})
    (root/'phase1-ketone-stereo-net.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)


if __name__=='__main__':
    run()
