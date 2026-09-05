import hashlib
import json
from collections import Counter
from pathlib import Path
import pytest
from cannabis_carbon.genome import _fasta
from cannabis_carbon.phase1_dopa_site_review import map_sites


def test_site_mapping_gaps_outside_and_mismatch():
    fields = [''] * 10 + ['2', '4', '2', '4', 'B-CD', 'BE-D']
    sites = map_sites(fields, 'ABCD', 'ABED', {1, 2, 3, 4})
    assert [s['status'] for s in sites] == ['outside-local-alignment', 'aligned-residue', 'query-gap', 'aligned-residue']
    assert sites[-1]['query_position'] == 4
    with pytest.raises(ValueError):
        map_sites(fields, 'AXCD', 'ABED', {2})


def test_pinned_dopa_sites_preserve_full_screen_and_feature_evidence(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    report = json.loads(Path('data/reports/phase1-dopa-site-review.json').read_text())
    search = json.loads(Path('data/reports/phase1-dopa-lyase-search.json').read_text())
    annotation = json.loads(Path('data/raw/phase1-dopa-lyase-search/Q3IWB0.json').read_text())
    assert report['model_eligible'] is False
    assert report['summary']['new_exact_enzyme_assignments'] == 0
    for p, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha
    expected_features = [f for f in annotation['features'] if f['type'] in ('Binding site', 'Active site', 'Mutagenesis', 'Modified residue')]
    assert report['reference_features'] == expected_features
    positions = {p for f in expected_features for p in range(f['location']['start']['value'], f['location']['end']['value'] + 1)}
    command = report['coordinate_command']
    lines = [l.split('\t') for l in Path(command[command.index('--out')+1]).read_text().splitlines()]
    assert len(lines) == report['coordinate_replay_alignment_count'] == 9
    assert Counter('\t'.join(f[:10]) for f in lines) == Counter(Path(search['hits_path']).read_text().splitlines())
    sequences = _fasta(Path(search['proteome_path']))
    assert len(sequences) == 30304
    assert {r['alignment_id'] for r in report['rows']} == {a['id'] for a in search['passing_alignments']}
    for row in report['rows']:
        f = row['alignment_columns']
        assert f in lines and row['model_eligible'] is False
        assert row['sites'] == map_sites(f, sequences[row['accession']], annotation['sequence']['value'], positions)
    assert report['summary']['reference_H89_aligned_residues'] == {'F': 7}
