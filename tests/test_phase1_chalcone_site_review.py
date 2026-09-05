import hashlib
import json
from collections import Counter
from pathlib import Path
from cannabis_carbon.genome import _fasta
from cannabis_carbon.phase1_dopa_site_review import map_sites


def test_chalcone_site_replay_preserves_all_hits_and_unmapped_leads(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-chalcone-site-review.json').read_text())
    search = json.loads(Path('data/reports/phase1-chalcone-search.json').read_text())
    annotation = json.loads(Path('data/raw/chalcone-annotations/P28012.json').read_text())
    assert report['model_eligible'] is False
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    assert report['reference_features'] == [f for f in annotation['features'] if f['type'] in ('Binding site', 'Active site', 'Site', 'Mutagenesis')]
    command = report['coordinate_command']
    lines = [l.split('\t') for l in Path(command[command.index('--out') + 1]).read_text().splitlines()]
    assert len(lines) == report['coordinate_replay_alignment_count'] == 115
    assert Counter('\t'.join(f[:10]) for f in lines) == Counter(Path(search['hits_path']).read_text().splitlines())
    queries = _fasta(Path(search['proteome_path']))
    refs = {r['accession']: r['sequence'] for r in search['reference_sequences']}
    assert len(queries) == 30304 and len(refs) == 50
    for f in lines:
        map_sites(f, queries[f[0].split('|')[1]], refs[f[1]], set())
    assert [r['alignment_columns'] for r in report['rows']] == [f for f in lines if f[1] == 'P28012']
    assert set(report['other_reference_leads_without_P28012_alignment']) == {r['accession'] for r in search['cannabis_candidates']} - {r['accession'] for r in report['rows']}
    for row in report['rows']:
        assert row['model_eligible'] is False
        assert row['sites'] == map_sites(row['alignment_columns'], queries[row['accession']], refs['P28012'], {48, 106, 113, 190})
        hits = [a for a in search['passing_alignments'] if a['cannabis_accession'] == row['accession'] and a['reference_accession'] == 'P28012']
        assert row['passes_original_screen'] == bool(hits)
        assert row['passing_alignment_id'] == (hits[0]['id'] if hits else None)
    passing = [r for r in report['rows'] if r['passes_original_screen']]
    assert len(passing) == 1
    assert passing[0]['accession'] == 'A0A7J6I409'
    assert [(s['query_position'], s['query_residue']) for s in passing[0]['sites']] == [(131, 'T'), (189, 'Y'), (196, 'N'), (273, 'S')]
