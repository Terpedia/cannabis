"""Validate source-protein annotations for ten exact MARTS stoichiometry gaps."""
import hashlib
import json
import re
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from rdkit import Chem, RDLogger
from .balance import _reaction_smiles_balance
from .phase1_balance_reference import carbon_participants
from .phase1_reference_discovery import direction_families

BOUNDARY = ('Source-protein annotation review, not Cannabis activity or corrected stoichiometry. '
    'Reaction references must match exact carbon-containing participants before evidence can be joined. '
    'Family, title and name-only matches do not resolve product identity. Atom tracing remains deferred.')


def queue(audit):
    targets = [t for t in audit['targets'] if t['marts_unbalanced_matches'] and not t['baseline_exact_balanced_participation']]
    ids = {m['source_record_id'] for t in targets for m in t['marts_unbalanced_matches']}
    return targets, [s for s in audit['source_ledger'] if s['id'] in ids]


def retrieve(accession, folder):
    path = folder / (accession + '.json')
    if path.exists():
        saved = json.loads(path.read_text())
        if saved['requested_accession'] != accession:
            raise ValueError('Cached accession mismatch')
        return saved
    result = {'requested_accession': accession}
    if accession.startswith('UPI'):
        result['status'] = 'UniParc-identifier-requires-separate-resolution'
        return result
    if not re.fullmatch(r'[A-Z0-9]{6,10}', accession):
        result['status'] = 'unsupported-accession-format'
        return result
    url = 'https://rest.uniprot.org/uniprotkb/' + accession + '.json'
    result['url'] = url
    try:
        with urllib.request.urlopen(url, timeout=45) as response:
            raw = response.read()
        entry = json.loads(raw)
        if accession not in [entry.get('primaryAccession')] + entry.get('secondaryAccessions', []):
            raise ValueError('Returned entry does not explicitly resolve requested accession')
        # Retain exact response bytes and their digest, plus parsed data for consumers.
        result.update(response_text=raw.decode(), response_sha256=hashlib.sha256(raw).hexdigest(),
            retrieved_at=datetime.now(timezone.utc).isoformat(), entry=entry)
        if entry.get('entryType') == 'Inactive':
            result['status'] = 'inactive-UniProt-entry-requires-explicit-resolution'
        else:
            sequence = entry.get('sequence', {}).get('value')
            if not sequence or len(sequence) != entry['sequence']['length']:
                raise ValueError('Missing or inconsistent sequence')
            result.update(status='retrieved', sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest())
    except (OSError, ValueError) as error:
        result.update(status='retrieval-failed', reason=str(error))
    if result['status'] in ('retrieved', 'inactive-UniProt-entry-requires-explicit-resolution'):
        path.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    return result


def stereo_relaxed(parts):
    """Diagnostic only: retain charge/isotopes/bonds/stoichiometry; never an evidence join."""
    if parts is None:
        return None
    result = []
    for side in parts:
        counts = Counter()
        for smiles, coefficient in side.items():
            mol = Chem.MolFromSmiles(smiles)
            Chem.RemoveStereochemistry(mol)
            counts[Chem.MolToSmiles(mol, isomericSmiles=True)] += coefficient
        result.append(dict(counts))
    return result


def build(audit, references, catalog, families):
    targets, sources = queue(audit)
    refs = {r['requested_accession']: r for r in references}
    if len(refs) != len(references):
        raise ValueError('Duplicate reference accession')
    for ref in references:
        if 'response_text' not in ref:
            if ref['status'] in ('retrieved', 'inactive-UniProt-entry-requires-explicit-resolution'):
                raise ValueError('Missing reference response provenance')
            continue
        raw = ref['response_text'].encode()
        if hashlib.sha256(raw).hexdigest() != ref['response_sha256'] or json.loads(raw) != ref['entry']:
            raise ValueError('Reference response checksum or parsed content mismatch')
        entry = ref['entry']
        if ref['requested_accession'] not in [entry.get('primaryAccession')] + entry.get('secondaryAccessions', []):
            raise ValueError('Unresolved accession substitution')
        if ref['status'] == 'retrieved':
            seq = entry['sequence']['value']
            if len(seq) != entry['sequence']['length'] or hashlib.sha256(seq.encode()).hexdigest() != ref['sequence_sha256']:
                raise ValueError('Reference sequence checksum mismatch')
    catalog_index = {r['rule_id']: r for r in catalog}
    if len(catalog_index) != len(catalog):
        raise ValueError('Ambiguous Rhea source identifiers')
    rows, retained_rhea, retained_families = [], {}, {}
    for source in sources:
        accession = source['source_record'].get('source_uniprot_id')
        ref = refs.get(accession)
        if accession and ref is None:
            raise ValueError('Missing reference retrieval outcome')
        entry = ref['entry'] if ref and ref['status'] == 'retrieved' else {}
        activities = [c for c in entry.get('comments', []) if c.get('commentType') == 'CATALYTIC ACTIVITY']
        matches = []
        original = carbon_participants(source['source_record']['reaction_smarts'])
        for i, activity in enumerate(activities):
            ids = {x['id'] for x in activity.get('reaction', {}).get('reactionCrossReferences', []) if x['database'] == 'Rhea'}
            for rid in sorted(ids):
                family = families.get(rid, {})
                if family:
                    retained_families[rid] = family
                candidates = sorted(({rid} | set(family.values())) & catalog_index.keys())
                assessments = []
                for candidate in candidates:
                    row = catalog_index[candidate]; retained_rhea[candidate] = row
                    parts = carbon_participants(row['reaction_smarts'])
                    element, charge = _reaction_smiles_balance(row['reaction_smarts'])
                    balanced = bool(element and charge and element['status'] == charge['status'] == 'balanced')
                    orientation = 'same-as-MARTS' if original and parts == original else 'reverse-of-MARTS' if original and parts == original[::-1] else None
                    relaxed_original, relaxed_parts = stereo_relaxed(original), stereo_relaxed(parts)
                    stereo_lead = bool(not orientation and relaxed_original and relaxed_parts in (relaxed_original, relaxed_original[::-1]))
                    assessments.append({'source_rhea_id': candidate, 'balanced': balanced,
                        'carbon_participant_match_orientation': orientation,
                        'stereo_only_diagnostic_lead': stereo_lead,
                        'diagnostic_boundary': 'Stereo-relaxed equality is not exact identity or enzyme evidence. Charge, isotopes, bond orders and coefficients remain unchanged.',
                        'status': 'exact-carbon-participants-balanced-alternative' if orientation and balanced else
                                  'exact-carbon-participants-balance-unresolved' if orientation else 'different-or-unresolved-carbon-participants'})
                matches.append({'activity_index': i, 'annotated_rhea_id': rid,
                    'family_mapping_status': 'explicit-published-family' if family else 'no-family-mapping',
                    'equation_assessments': assessments})
        exact = any(a['status'] == 'exact-carbon-participants-balanced-alternative' for m in matches for a in m['equation_assessments'])
        rows.append({'id': source['id'], 'source_record': source,
            'target_ids': [t['cannabisdb_id'] for t in targets if any(m['source_record_id'] == source['id'] for m in t['marts_unbalanced_matches'])],
            'requested_accession': accession, 'source_genbank_id': source['source_record'].get('source_genbank_id'),
            'reference_status': ref['status'] if ref else 'no-UniProt-reference-in-source',
            'resolved_accession': entry.get('primaryAccession'), 'source_organism': entry.get('organism'),
            'source_entry_type': entry.get('entryType'), 'catalytic_activity_reviews': matches,
            'status': 'exact-balanced-reference-alternative-found-review-required' if exact else
                'source-annotations-do-not-resolve-exact-equation' if activities else 'source-catalytic-annotation-unresolved',
            'enzyme_evidence_ids': [], 'claim_boundary': BOUNDARY})
    return {'schema': 'cannabis-carbon.phase1-marts-gap-references.v1', 'targets': targets, 'rows': rows,
        'references': references, 'rhea_source_records': list(retained_rhea.values()),
        'rhea_direction_families': retained_families, 'claim_boundary': BOUNDARY,
        'summary': {'target_records': len(targets), 'MARTS_source_rows': len(rows), 'reference_identifiers': len(references),
            'reference_status_counts': dict(Counter(r['status'] for r in references)),
            'source_status_counts': dict(Counter(r['status'] for r in rows)),
            'sources_with_catalytic_Rhea_annotations': sum(bool(r['catalytic_activity_reviews']) for r in rows),
            'sources_with_stereo_only_diagnostic_leads': sum(any(a['stereo_only_diagnostic_lead'] for m in r['catalytic_activity_reviews'] for a in m['equation_assessments']) for r in rows),
            'exact_balanced_reference_alternatives': sum(r['status'] == 'exact-balanced-reference-alternative-found-review-required' for r in rows)},
        'next_tests': ['Review cited publications and supplements for each exact substrate/product, including minor products and stereochemistry.',
            'Resolve UniParc and GenBank-only references explicitly; preserve sequence and identifier provenance.',
            'Curate complete stoichiometry before transferring reference evidence or screening exact-reaction candidates.']}


def run():
    RDLogger.DisableLog('rdApp.warning'); RDLogger.DisableLog('rdApp.error')
    audit_path = Path('data/reports/phase1-marts-audit.json')
    catalog_path = Path('data/raw/phase1-balance-reference-catalog.json')
    family_path = Path('data/raw/phase1-reference-discovery/rhea-directions.tsv')
    audit = json.loads(audit_path.read_text())
    _, sources = queue(audit)
    accessions = sorted({s['source_record']['source_uniprot_id'] for s in sources if s['source_record'].get('source_uniprot_id')})
    folder = Path('data/raw/phase1-marts-gap-references'); folder.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        references = list(pool.map(lambda a: retrieve(a, folder), accessions))
    report = build(audit, references, json.loads(catalog_path.read_text()), direction_families(family_path.read_text()))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (audit_path, catalog_path, family_path)}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-marts-gap-references.json').write_text(payload)
    Path('docs/data/phase1-marts-gap-references.json').write_text(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('target', report['targets']), ('source_review', report['rows']), ('reference', references),
              ('rhea_source', report['rhea_source_records'])]
    metadata = {k: v for k, v in report.items() if k not in ('targets', 'rows', 'references', 'rhea_source_records')}
    with Path('data/derived/phase1-marts-gap-references.ndjson').open('w') as handle:
        for kind, records in [('metadata', [metadata])] + groups:
            for record in records:
                rid = record.get('id', record.get('cannabisdb_id', record.get('requested_accession', record.get('rule_id', 'metadata'))))
                handle.write(json.dumps({'record_kind': kind, 'record_id': rid,
                    'record_json': json.dumps(record, separators=(',', ':')), 'report_sha256': digest}) + '\n')
    print(json.dumps({'sha256': digest, 'bytes': len(payload.encode()), **report['summary']}))


if __name__ == '__main__':
    run()
