import copy
import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_completion_protein_discovery import build
from cannabis_carbon.phase1_completion_protein_evidence import build as evidence_build


def test_exact_equation_join_rejects_changed_coefficients():
    c = {'completions': [{'id': 'h', 'balanced_equation_id': 'r', 'left': [{'compound_id': 'a', 'coefficient': 1}], 'right': []}], 'variants': [], 'targets': []}
    with pytest.raises(ValueError, match='stoichiometry'):
        build(c, {'source_ledger': []}, {'candidate_reaction_evidence_ids': {'r': ['e']}},
              {'reactions': [{'id': 'r', 'left': [], 'right': []}]}, {'references': [], 'rows': []})


def test_published_discovery_replays_original_sources_and_keeps_exclusions():
    root = Path(__file__).resolve().parents[1]
    report = json.loads((root / 'data/reports/phase1-completion-protein-discovery.json').read_text())
    inputs = []
    for name, digest in report['source_sha256'].items():
        raw = (root / name).read_bytes(); assert hashlib.sha256(raw).hexdigest() == digest
        inputs.append(json.loads(raw))
    before = copy.deepcopy(inputs)
    assert build(*inputs) == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert inputs == before
    assert len(report['rows']) + len(report['existing_evidence_matches']) == 765
    for row in report['rows']:
        assert all(set(ref['source_record_ids']) <= set(row['source_record_ids']) for ref in row['reference_matches'])


def test_published_overlay_preserves_chemistry_and_uncertainty():
    root = Path(__file__).resolve().parents[1]
    raw = (root / 'data/reports/phase1-completion-protein-evidence.json').read_bytes()
    assert raw == (root / 'docs/data/phase1-completion-protein-evidence.json').read_bytes()
    report = json.loads(raw); inputs = []
    for name, digest in report['source_sha256'].items():
        data = (root / name).read_bytes(); assert hashlib.sha256(data).hexdigest() == digest
        inputs.append(json.loads(data))
    before = copy.deepcopy(inputs)
    assert evidence_build(*inputs) == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert inputs == before
    for row in inputs[2]['rows']:
        assert 'inferred-inorganic-stoichiometry-unverified' in row['validation_blockers']
        assert 'original-MARTS-exact-product-identity-unverified' in row['validation_blockers']
        if row['screened_cannabis_proteins']:
            assert row['evidence_class'] == 'original-MARTS-source-homology-for-inferred-stoichiometry'
