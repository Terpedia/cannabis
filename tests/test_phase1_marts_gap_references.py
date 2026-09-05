import copy
import hashlib
import json
from pathlib import Path
import pytest
from cannabis_carbon.phase1_marts_gap_references import build, stereo_relaxed
from cannabis_carbon.phase1_balance_reference import carbon_participants


def reference():
    entry = {'primaryAccession': 'P12345', 'sequence': {'value': 'ACDE', 'length': 4},
        'comments': [{'commentType': 'CATALYTIC ACTIVITY', 'reaction': {'reactionCrossReferences': [
            {'database': 'Rhea', 'id': 'RHEA:101'}]}}]}
    raw = json.dumps(entry)
    return {'requested_accession': 'P12345', 'status': 'retrieved', 'entry': entry,
        'response_text': raw, 'response_sha256': hashlib.sha256(raw.encode()).hexdigest(),
        'sequence_sha256': hashlib.sha256(b'ACDE').hexdigest()}


def audit(smiles='C=C>>CCO'):
    return {'targets': [{'cannabisdb_id': 'T', 'baseline_exact_balanced_participation': False,
                        'marts_unbalanced_matches': [{'source_record_id': 'S'}]}],
        'source_ledger': [{'id': 'S', 'source_record': {'source_uniprot_id': 'P12345', 'reaction_smarts': smiles}}]}


def test_explicit_family_and_exact_participants_required_for_reference_join():
    # Deliberately nonadjacent ID: no arithmetic family inference is permitted.
    families = {'RHEA:101': {'RHEA_ID_MASTER': 'RHEA:101', 'RHEA_ID_LR': 'RHEA:777'}}
    catalog = [{'rule_id': 'RHEA:777', 'reaction_smarts': 'C=C.O>>CCO'}]
    result = build(audit(), [reference()], catalog, families)
    assert result['summary']['exact_balanced_reference_alternatives'] == 1
    assert result['rows'][0]['enzyme_evidence_ids'] == []
    assert build(audit(), [reference()], catalog, {})['summary']['exact_balanced_reference_alternatives'] == 0
    changed = [{'rule_id': 'RHEA:777', 'reaction_smarts': 'C=C>>CCO'}]
    assert build(audit(), [reference()], changed, families)['summary']['exact_balanced_reference_alternatives'] == 0


def test_stereo_relaxation_is_diagnostic_and_retains_charge_and_isotopes():
    left = 'C[C@H](O)Cl'; right = 'C[C@@H](O)Cl'
    result = build(audit(left + '>>' + left), [reference()],
        [{'rule_id': 'RHEA:101', 'reaction_smarts': right + '>>' + right}], {})
    assert result['summary']['sources_with_stereo_only_diagnostic_leads'] == 1
    assert result['summary']['exact_balanced_reference_alternatives'] == 0
    for a, b in [('CO>>CO', 'C[O-]>>C[O-]'), ('[13CH3]O>>CO', 'CO>>CO')]:
        assert stereo_relaxed(carbon_participants(a)) != stereo_relaxed(carbon_participants(b))


def test_reference_provenance_fails_closed():
    for key, value in [('response_sha256', 'bad'), ('sequence_sha256', 'bad'), ('requested_accession', 'Q12345')]:
        ref = reference(); ref[key] = value
        with pytest.raises(ValueError):
            build(audit(), [ref], [], {})
    ref = reference(); del ref['response_text']
    with pytest.raises(ValueError, match='provenance'):
        build(audit(), [ref], [], {})


def test_published_reference_review_retains_sources_and_replays_without_raw_cache():
    root = Path(__file__).resolve().parents[1]
    raw = (root / 'data/reports/phase1-marts-gap-references.json').read_bytes()
    report = json.loads(raw)
    assert raw == (root / 'docs/data/phase1-marts-gap-references.json').read_bytes()
    for name, digest in report['source_sha256'].items():
        path = root / name
        if path.exists():
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    source = json.loads((root / 'data/reports/phase1-marts-audit.json').read_text())
    before = copy.deepcopy(source)
    rebuilt = build(source, report['references'], report['rhea_source_records'], report['rhea_direction_families'])
    assert rebuilt == {k: v for k, v in report.items() if k != 'source_sha256'}
    assert source == before
    assert len(report['targets']) == 10 and len(report['rows']) == 32
    assert len(report['references']) == 27
    assert all(r['enzyme_evidence_ids'] == [] for r in report['rows'])
