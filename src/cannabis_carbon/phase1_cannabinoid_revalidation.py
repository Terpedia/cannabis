"""Revalidate legacy synthase reference leads against the pinned whole proteome."""
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from .genome import _fasta
from .phase1_new_protein_search import run as search, export_table

SOURCE_IDS = ('cannabis:reaction:cbga-to-cbda', 'cannabis:reaction:cbga-to-thca')


def prepare(audit, sequences, references):
    rows, checks = [], []
    for source_id in SOURCE_IDS:
        matches = [(r, s) for r in audit['rows'] for s in r['core_source_reactions'] if s['id'] == source_id]
        if len(matches) != 1:
            raise ValueError('Ambiguous core equation link')
        equation, source = matches[0]
        refs = []
        for protein in source['candidate_proteins']:
            old = protein.get('specialized_search') or {}
            accession = old.get('reference_accession')
            if not accession:
                continue
            record = references[accession]
            if record['primaryAccession'] != accession:
                raise ValueError('Reference accession redirected')
            ec = {e['value'] for e in record['proteinDescription']['recommendedName'].get('ecNumbers', [])}
            if not set(source['ec_numbers']) <= ec:
                raise ValueError('Legacy source EC not corroborated by reference record')
            seq = sequences[protein['accession']]
            checks.append({'source_reaction_id': source_id, 'candidate_accession': protein['accession'],
                'pinned_length': len(seq), 'pinned_sequence_sha256': hashlib.sha256(seq.encode()).hexdigest(),
                'legacy_specialized_search': old,
                'legacy_query_lengths_match_pinned_sequence': all(h['query_length'] == len(seq) for h in old['hits']),
                'action': 'Retain legacy assertion; use new whole-proteome search for current alignment evidence.'})
            refs.append({'accession': accession, 'join_method': 'legacy-core-reference-lead-with-current-exact-EC-corroboration',
                'ec_numbers': sorted(ec), 'reference_record': record,
                'claim_boundary': 'Reference annotation and legacy source linkage, not validation of every encoded substrate/product stereoisomer or the candidate protein activity.'})
        if not refs:
            raise ValueError('No linked synthase reference')
        rows.append({'reaction_id': equation['reaction_id'], 'left': equation['reaction']['left'],
            'right': equation['reaction']['right'], 'sources': equation['reaction']['sources'],
            'target_ids': equation['target_ids'], 'priority_target_ids': equation['target_ids'],
            'hypothesis_ids': [], 'reference_matches': refs,
            'legacy_source_reaction_id': source_id})
    return {'schema': 'cannabis-carbon.phase1-cannabinoid-revalidation-references.v1', 'rows': rows,
        'legacy_checks': checks,
        'claim_boundary': 'Revalidation of two specific legacy reference links. Exact reaction structures and source orientation remain unpromoted. Sequence homology is not assay evidence or a complete CO2 pathway.'}


def run():
    raw = Path('data/raw/phase1-cannabinoid-revalidation'); raw.mkdir(parents=True, exist_ok=True)
    audit_path = Path('data/reports/phase1-producer-screen-audit.json')
    audit = json.loads(audit_path.read_text())
    for path, sha in audit['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
            raise ValueError('Audit lineage changed')
    accessions = sorted({p['specialized_search']['reference_accession'] for r in audit['rows']
        for s in r['core_source_reactions'] if s['id'] in SOURCE_IDS for p in s['candidate_proteins']})
    references, lookups = {}, []
    for accession in accessions:
        url = f'https://rest.uniprot.org/uniprotkb/{accession}.json'
        path = raw / (accession + '.json')
        receipt = raw / (accession + '.request.json')
        if not receipt.exists():
            with urllib.request.urlopen(url, timeout=45) as response:
                payload = response.read()
            json.loads(payload)
            path.write_bytes(payload)
            receipt.write_text(json.dumps({'url': url, 'retrieved_at': datetime.now(timezone.utc).isoformat(),
                'snapshot': str(path), 'sha256': hashlib.sha256(payload).hexdigest()}) + '\n')
        lookup = json.loads(receipt.read_text())
        if lookup['url'] != url or hashlib.sha256(path.read_bytes()).hexdigest() != lookup['sha256']:
            raise ValueError('Reference snapshot mismatch')
        references[accession] = json.loads(path.read_text()); lookups.append(lookup)
    proteome = Path('data/raw/UP000583929.fasta')
    report = prepare(audit, _fasta(proteome), references)
    report['lookups'] = lookups
    paths = [audit_path, proteome, *(Path(l['snapshot']) for l in lookups)]
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    discovery = Path('data/reports/phase1-cannabinoid-revalidation-references.json')
    discovery.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    output = Path('data/reports/phase1-cannabinoid-revalidation-search.json')
    search(discovery, raw, output, evidence_class='legacy-reference-revalidated-homology; exact-chemistry-unverified',
        additional_blockers=('legacy-encoded-substrate-product-specificity-unverified',))
    export_table(output, Path('data/derived/phase1-cannabinoid-revalidation-search.ndjson'))


if __name__ == '__main__':
    run()
