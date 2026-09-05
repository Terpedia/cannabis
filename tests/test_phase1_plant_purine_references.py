import copy
import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_plant_purine_references import build
from cannabis_carbon.phase1_reference_discovery import direction_families

ROOT=Path(__file__).resolve().parents[1]


def test_plant_reference_audit_replays_without_promoting_unreviewed_annotations():
    def read(n):return json.loads((ROOT/'data/reports'/(n+'.json')).read_text())
    network,expanded,report,search=[read(n) for n in ('phase1-full-balanced-network','phase1-expanded-candidate-net','phase1-plant-purine-references','phase1-plant-purine-search')]
    before=copy.deepcopy((network,expanded))
    families=direction_families((ROOT/'data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text())
    tsv=(ROOT/report['lookup']['snapshot']).read_text()
    assert build(tsv,network,expanded,families)=={k:v for k,v in report.items() if k not in ('lookup','source_sha256')}
    assert (network,expanded)==before
    for path,sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha
    assert report['summary']['review_status_counts']=={'reviewed':11,'unreviewed':13}
    assert report['summary']['family_status_counts']=={'already-candidate-linked':10,'candidate-evidence-gap':3}
    assert report['summary']['references_without_Rhea_annotation']==['Q9SJ42']
    assert {r['source_master_id'] for r in report['rows']}=={'RHEA:16853','RHEA:23920','RHEA:22192'}
    assert all(p['review_status']=='unreviewed' for r in report['rows'] for p in r['reference_matches'])
    assert all(r['reaction_id'] not in expanded['candidate_reaction_evidence_ids'] for r in report['rows'])
    assert search['summary']['proteome_sequences']==30304
    assert search['summary']['retrieved_references']==8
    assert search['summary']['distinct_cannabis_candidates']==4
    assert search['summary']['equations_with_screened_candidates']==3
    assert search['summary']['passing_alignments']==16
    assert all('unreviewed' in r['evidence_class'] and 'unreviewed-reference-annotation' in r['validation_blockers'] for r in search['rows'])
