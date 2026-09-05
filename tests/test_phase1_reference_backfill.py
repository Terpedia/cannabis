import hashlib
import json
from pathlib import Path

from cannabis_carbon.phase1_reference_backfill import queue
from cannabis_carbon.phase1_reference_discovery import direction_families
from cannabis_carbon.phase1_new_references import attach

ROOT = Path(__file__).resolve().parents[1]


def test_backfill_replays_prior_lookup_coverage_and_exact_reference_joins(monkeypatch):
    monkeypatch.chdir(ROOT)
    read = lambda p: json.loads(Path(p).read_text())
    report = read('data/reports/phase1-reference-backfill.json')
    for p, digest in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == digest
    prior_paths = ['data/reports/'+n+'.json' for n in ('phase1-new-references','phase1-route-references','phase1-catalog-references')]
    lookups = [l for p in prior_paths for l in read(p)['lookups']]
    screens = {p:read(p) for p in report['source_sha256'] if p.endswith('protein-search.json')}
    families = direction_families(Path(report['rhea_direction_source']['snapshot']).read_text())
    rows, audit = queue(read('data/reports/phase1-catalog-net-gaps.json'),
        read('data/reports/phase1-catalog-evidence.json'), screens, lookups, families)
    assert audit == report['remaining_gap_audit']
    assert len(audit) == 368 and len(rows) == 157
    assert all(not p['reference_sequences_present'] for a in audit
        if a['id'] in {r['reaction_id'] for r in rows} for p in a['prior_screens'])
    proteins = attach(rows,report['lookups'])
    for row in rows:
        if not row['rhea_families']:
            row['lookup_status'] = 'no-published-family-mapping'
    assert rows == report['rows']
    used = {m['accession'] for row in rows for m in row['reference_matches']}
    assert [proteins[a] for a in sorted(used)] == report['reference_proteins']
    previous = {m for l in lookups if l['status']=='retrieved' for m in l['requested_master_ids']}
    new = {m for l in report['lookups'] if 'reused_from' not in l for m in l['requested_master_ids']}
    assert len(new) == 51 and not new & previous


def test_not_queried_is_distinct_from_completed_no_reference_lookup():
    gap = {'id':'r','reaction_id':'r','selected_certificate_target_ids':['t']}
    reaction = {'id':'r','left':[],'right':[],'sources':[{'source_reaction_id':'RHEA:1'}]}
    catalog = {'gap_priorities':[gap],'reactions':[reaction]}
    row = {'reaction_id':'r','passing_alignment_ids':[],'search_status':'no-reference-sequence',
        'reference_sequences_present':[],'reference_sequences_missing':[],'raw_alignment_count':0}
    screens = {'prior':{'rows':[row]}}
    families = {'RHEA:1':{'RHEA_ID_MASTER':'RHEA:10'}}
    args = (catalog,{'enzyme_evidence':[]},screens)
    assert queue(*args,[],families)[1][0]['reference_gap_status']=='reference-family-not-queried'
    completed = [{'status':'retrieved','requested_master_ids':['RHEA:10']}]
    assert queue(*args,completed,families)[1][0]['reference_gap_status']=='reviewed-family-lookup-completed'
    assert queue(*args,completed,{})[1][0]['reference_gap_status']=='no-published-family-mapping'
    row['reference_sequences_present']=['protein']
    assert queue(*args,completed,families)[0]==[]


def test_backfill_candidates_reduce_reaction_gaps_without_closing_more_certificates():
    read = lambda n: json.loads((ROOT / 'data/reports' / (n+'.json')).read_text())
    catalog, previous, search = [read(n) for n in ('phase1-catalog-net-gaps','phase1-catalog-evidence','phase1-backfill-protein-search')]
    old = {e['reaction_id'] for e in previous['enzyme_evidence']}
    new = {r['reaction_id'] for r in search['rows'] if r['screened_cannabis_proteins']}
    assert len(new)==19 and not old & new
    assert len({g['id'] for g in catalog['gap_priorities']} - old - new)==349
    additionally_closed = [t['cannabisdb_id'] for t in catalog['targets']
        if t['missing_candidate_reaction_ids'] and not set(t['missing_candidate_reaction_ids']) <= old
        and set(t['missing_candidate_reaction_ids']) <= old | new]
    assert additionally_closed == []
