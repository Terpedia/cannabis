import json
import pytest
from cannabis_carbon.phase1_overlay import build


def test_exact_variant_scope_and_candidate_not_confirmation(tmp_path):
    queue, evidence, output = [tmp_path / n for n in ('queue', 'evidence', 'output')]
    rows = [{'reaction_id': 'RHEA:1', 'reaction_smarts': s, 'balance_status': b, 'source_urls': []}
            for s, b in [('a>>b', 'balanced'), ('c>>d', 'balanced'), ('e>>f', 'imbalanced')]]
    queue.write_text(json.dumps({'rows': rows}))
    evidence.write_text(json.dumps({'rows': [
        {**rows[0], 'search_status': 'hits-found', 'sequence_hits': [{'passes_screen': True, 'cannabis_accession': 'P1'}]},
        {**rows[1], 'search_status': 'no-reference-sequence', 'core_reaction_evidence': [
            {'core_reaction_id': 'rhea:1', 'enzyme_association_ids': ['P2']}]}]}))
    summary = build(queue, evidence, output)
    assert summary['balanced_variants'] == 2
    assert summary['balanced_variants_with_candidate_enzyme_evidence'] == 2
    assert summary['evidence_status_counts'] == {'screened-homology': 1, 'core-association': 1, 'balance-unresolved': 1}
    assert json.loads(output.read_text())['rows'][2]['search_status'].startswith('not-searched')
    evidence.write_text(json.dumps({'rows': []}))
    with pytest.raises(ValueError, match='exactly'):
        build(queue, evidence, output)


def test_duplicate_variants_rejected(tmp_path):
    queue, evidence, output = [tmp_path / n for n in ('queue', 'evidence', 'output')]
    row = {'reaction_id': 'RHEA:1', 'reaction_smarts': 'a>>b'}
    queue.write_text(json.dumps({'rows': [row, row]}))
    evidence.write_text(json.dumps({'rows': []}))
    with pytest.raises(ValueError, match='Duplicate'):
        build(queue, evidence, output)
