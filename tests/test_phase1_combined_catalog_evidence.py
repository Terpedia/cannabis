import copy
import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_combined_catalog_evidence import build

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT/'data/reports'/(name+'.json')).read_text())


def test_combined_evidence_replays_exact_union_and_frozen_chemistry():
    catalog,previous,search,report = [read(n) for n in ('phase1-catalog-net-gaps','phase1-catalog-evidence',
        'phase1-backfill-protein-search','phase1-combined-catalog-evidence')]
    before = copy.deepcopy((catalog,previous,search))
    assert build(catalog,previous,search)=={k:v for k,v in report.items() if k!='source_sha256'}
    assert (catalog,previous,search)==before
    for path,sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha
    old={e['reaction_id']:e for e in previous['enzyme_evidence']}
    combined={e['reaction_id']:e for e in report['enzyme_evidence']}
    assert len(combined)==116
    assert all(combined[rid]==e for rid,e in old.items())
    expected={r['reaction_id'] for r in search['rows'] if r['screened_cannabis_proteins']}
    assert combined.keys()-old.keys()==expected
    for rid in expected:
        assert combined[rid]['full_search_report'].endswith('phase1-backfill-protein-search.json')
        assert combined[rid]['evidence_class']=='reference-backfill-direction-unresolved-homology-candidate'
    updates={c['compound_id']:c for c in report['certificate_updates']}
    for cert in catalog['certificates']:
        before_missing=cert['missing_candidate_reaction_ids']
        remaining=[rid for rid in before_missing if rid not in combined]
        if remaining!=before_missing:
            assert updates[cert['compound_id']]['missing_candidate_reaction_ids']==remaining
            assert updates[cert['compound_id']]['baseline_missing_candidate_reaction_ids']==before_missing
        else:
            assert cert['compound_id'] not in updates
    assert report['summary']['remaining_missing_candidate_equations']==349
    assert report['summary']['selected_certificate_targets_with_candidates_for_all_steps']==102
    assert report['summary']['selected_certificate_targets_with_remaining_gaps']==202
    payload=(ROOT/'data/reports/phase1-combined-catalog-evidence.json').read_bytes()
    assert payload==(ROOT/'docs/data/catalog-net-view/evidence.json').read_bytes()
    manifest=json.loads((ROOT/'docs/data/catalog-net-view/index.json').read_text())
    assert manifest['evidence']['sha256']==hashlib.sha256(payload).hexdigest()
    assert manifest['evidence']['bytes']==len(payload)


def test_combined_supplement_rejects_repeated_evidence():
    catalog,previous,search=[read(n) for n in ('phase1-catalog-net-gaps','phase1-catalog-evidence','phase1-backfill-protein-search')]
    previous['enzyme_evidence'].append(previous['enzyme_evidence'][0])
    with pytest.raises(ValueError,match='Duplicate'):
        build(catalog,previous,search)
