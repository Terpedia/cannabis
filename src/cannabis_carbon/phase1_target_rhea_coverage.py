"""Expand the whole-target participation audit to Terpedia's full Rhea snapshot."""
import hashlib
import json
from pathlib import Path
from rdkit import RDLogger, rdBase
from .phase1_target_coverage import audit_targets, export_table


def source_records(rows):
    return [{
        'id': 'terpedia-full-rhea:' + hashlib.sha256(json.dumps(
            [r['rule_id'], r['reaction_smarts']], separators=(',', ':')).encode()).hexdigest(),
        'source_reaction_id': r['rule_id'], 'reaction_smiles': r['reaction_smarts'],
        'source_url': r['source_url'], 'source_download_url': r['source_download_url'],
        'source_evidence_type': r['source_evidence_type'],
        'source_layer': 'terpedia-full-rhea-catalog',
        'direction_evidence': {'source_direction_mode': r['direction_mode']},
        'orientation_boundary': 'As-recorded source equation orientation, not established Cannabis physiological direction.',
        'organism_evidence': 'No Cannabis enzyme activity inferred from a Rhea structure match.'
    } for r in rows]


def run():
    RDLogger.DisableLog('rdApp.warning')
    RDLogger.DisableLog('rdApp.error')
    metadata_path = Path('data/reports/phase1-balance-reference.json')
    metadata = json.loads(metadata_path.read_text())['catalog']
    raw = Path(metadata['snapshot'])
    if hashlib.sha256(raw.read_bytes()).hexdigest() != metadata['sha256']:
        raise ValueError('Terpedia Rhea snapshot checksum mismatch')
    source = json.loads(raw.read_text())
    if len(source) != metadata['row_count']:
        raise ValueError('Terpedia Rhea snapshot row count mismatch')
    paths = [Path('docs/data/compounds.json'), Path('docs/data/networkdb.json'),
             Path('data/reports/phase1-reaction-catalog.json')]
    targets, network, catalog = [json.loads(p.read_text()) for p in paths]
    report = audit_targets(targets['compounds'], network, catalog,
                           source_records(source), matching_ledger_only=True)
    report['catalog_provenance'] = metadata
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in paths + [raw, metadata_path]}
    report['rdkit_version'] = rdBase.rdkitVersion
    report['summary']['retained_matching_reaction_records'] = len(report['reaction_ledger'])
    report['summary']['additional_rhea_source_records'] = len(source)
    report['metric_scope'] += ' Full Rhea includes reactions from other organisms. Overlapping source layers are not deduplicated reaction counts. Ledger retains matching equations only; balance counts cover every scanned source record.'
    baseline_path = Path('data/reports/phase1-target-coverage.json')
    baseline = json.loads(baseline_path.read_text())
    if [r['cannabisdb_id'] for r in baseline['targets']] != [r['cannabisdb_id'] for r in report['targets']]:
        raise ValueError('Baseline target inventory mismatch')
    report['baseline_comparison'] = {
        'path': str(baseline_path), 'sha256': hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
        'new_balanced_participants': [r['cannabisdb_id'] for old, r in zip(baseline['targets'], report['targets'])
                                      if old['coverage_status'] != 'balanced-reaction-participant'
                                      and r['coverage_status'] == 'balanced-reaction-participant']}
    output = Path('data/reports/phase1-target-rhea-coverage.json')
    output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    count = export_table(output, output.with_suffix('.ndjson'))
    print(json.dumps({'summary': report['summary'], 'new_balanced_participants':
                     len(report['baseline_comparison']['new_balanced_participants']), 'export_rows': count}))


if __name__ == '__main__':
    run()
