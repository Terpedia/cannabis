import csv
import json
from pathlib import Path
import pytest
from rdkit import Chem
from cannabis_carbon.phase1_chalcone_reference_search import build, RID
from cannabis_carbon.phase1_screened_overlay import build_overlay
from test_phase1_new_protein_search import test_published_search_preserves_gap_scope_thresholds_sequences_and_joins as verify_search


def test_chalcone_discovery_preserves_exact_stereo_and_reference_records(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-chalcone-references.json').read_text())
    assert report == build(report['lookups'][0])
    gap = next(r for r in json.loads(Path('data/reports/phase1-current-gap-priority.json').read_text())['rows'] if r['reaction_id'] == RID)
    row = report['rows'][0]
    for side in ('left', 'right'):
        assert row[side] == gap['reaction'][side]
    assert row['participants'] == gap['participants']
    product = next(p for p in row['participants'] if p['id'] == row['right'][0]['compound_id'])
    assert [s for _, s in Chem.FindMolChiralCenters(Chem.MolFromSmiles(product['smiles']))] == ['S']
    assert row['target_ids'] == gap['remaining_target_ids']
    assert len(row['target_ids']) == 7
    assert row['historical_selected_uses'] == gap['selected_uses']
    assert {t for u in row['selected_uses'] for t in u['target_ids']} == set(row['target_ids'])
    item = report['lookups'][0]
    with Path(item['snapshot']).open() as stream:
        records = list(csv.DictReader(stream, delimiter='\t'))
    assert len(records) == report['summary']['reference_leads'] == 50
    assert report['summary']['references_without_rhea_field'] == 50
    assert {m['accession']: m['source_records'][0]['record'] for m in row['reference_matches']} == {r['Entry']: r for r in records}
    assert all(m['model_eligible'] is False and m['exact_reaction_annotation_match'] is False for m in row['reference_matches'])


def test_chalcone_full_proteome_screen_retains_weak_hits_and_blocks_promotion(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    verify_search('phase1-chalcone-search')
    search = json.loads(Path('data/reports/phase1-chalcone-search.json').read_text())
    assert search['summary']['proteome_sequences'] == 30304
    assert search['summary']['retrieved_references'] == search['summary']['requested_references'] == 50
    assert search['summary']['raw_alignments'] == 115
    assert search['summary']['passing_alignments'] == 60
    assert search['summary']['distinct_cannabis_candidates'] == 5
    assert search['model_eligible'] is False
    assert all(r['model_eligible'] is False for r in search['rows'])
    with pytest.raises(ValueError, match='ineligible'):
        build_overlay({'reactions': [{'id': RID}], 'hypotheses': []}, search)
