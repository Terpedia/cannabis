"""Screen new exact-source Rhea-family reference leads against the full proteome."""
import hashlib
import json
import math
import subprocess
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from .genome import _fasta
from .phase1_search import fasta_accession
from .phase1_family_search import parse_references, parse_hits


def annotate(discovery, references, hits):
    passing = {}
    for alignments in hits.values():
        for hit in alignments:
            if hit['passes_screen']:
                aid = hashlib.sha256(json.dumps(hit, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
                passing[aid] = {'id': aid, **hit}
    by_reference = {}
    for aid, hit in passing.items():
        by_reference.setdefault(hit['reference_accession'], []).append(aid)
    rows = []
    for gap in discovery['rows']:
        requested = {r['accession'] for r in gap['reference_matches']}
        present = requested & references.keys()
        ids = sorted({aid for ref in present for aid in by_reference.get(ref, [])})
        raw_count = sum(len(hits.get(ref, [])) for ref in present)
        proteins = sorted({passing[aid]['cannabis_accession'] for aid in ids})
        rows.append({'reaction_id': gap['reaction_id'], 'target_ids': gap['target_ids'],
            'priority_target_ids': gap['priority_target_ids'], 'hypothesis_ids': gap['hypothesis_ids'],
            'reference_matches': gap['reference_matches'],
            'reference_sequences_present': sorted(present), 'reference_sequences_missing': sorted(requested - present),
            'passing_alignment_ids': ids, 'raw_alignment_count': raw_count,
            'screened_cannabis_proteins': proteins,
            'search_status': 'screened-candidates' if proteins else 'weak-hits-only' if raw_count else 'no-hits' if present else 'no-reference-sequence',
            'evidence_class': 'direction-unresolved-reference-homology-candidate' if proteins else 'unresolved',
            'validation_blockers': ['physiological-direction-unverified', 'exact-substrate-specificity-unverified',
                'catalytic-residues-and-domains-not-reviewed', 'compartment-and-expression-unverified', 'all-input-supply-unestablished'],
            'proposed_test': 'Review exact reference reactions, domains and catalytic residues; assay candidate proteins using every required input in the parent hypothesis, with authentic product standards, a characterized reference and no-enzyme controls. Establish direction and tissue compatibility before claiming a Cannabis pathway.'})
    return rows, list(passing.values())


def run():
    source = Path('data/reports/phase1-new-references.json')
    discovery = json.loads(source.read_text())
    for filename, digest in discovery['source_sha256'].items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest() != digest:
            raise ValueError('Discovery source checksum mismatch')
    ids = sorted({m['accession'] for r in discovery['rows'] for m in r['reference_matches']})
    raw = Path('data/raw/phase1-new-protein-search'); raw.mkdir(parents=True, exist_ok=True)
    def retrieve(batch):
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({
            'query': ' OR '.join('accession:' + accession for accession in batch), 'format': 'fasta'})
        retrieval = {'requested_accessions': batch, 'url': url}
        parsed = {}
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                data = response.read()
            parsed = parse_references(data, set(batch))
            digest = hashlib.sha256(data).hexdigest()
            snapshot = raw / (digest + '.fasta'); snapshot.write_bytes(data)
            retrieval.update(status='retrieved', sha256=digest, snapshot=str(snapshot), missing_accessions=sorted(set(batch) - parsed.keys()))
        except (OSError, ValueError) as error:
            retrieval.update(status='retrieval-failed', reason=str(error))
        return retrieval, parsed
    batches = [ids[start:start + 60] for start in range(0, len(ids), 60)]
    retrievals, references = [], {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        for n, (retrieval, parsed) in enumerate(pool.map(retrieve, batches), 1):
            retrievals.append(retrieval); references.update(parsed)
            print(f'Reference batch {n}/{len(batches)}: {retrieval["status"]}', flush=True)
    if not references:
        raise ValueError('No reference sequences retrieved')
    reference_path = raw / 'references.fasta'
    reference_path.write_text(''.join(f'>{acc} {r["header"]}\n{r["sequence"]}\n' for acc, r in sorted(references.items())))
    (raw / 'retrievals.json').write_text(json.dumps({'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'retrievals': retrievals}, separators=(',', ':')) + '\n')
    proteome = Path('data/raw/UP000583929.fasta')
    prior = json.loads(Path('data/reports/phase1-family-protein-search.json').read_text())
    if hashlib.sha256(proteome.read_bytes()).hexdigest() != prior['proteome_sha256']:
        raise ValueError('Cannabis proteome differs from verified prior snapshot')
    queries = _fasta(proteome)
    headers = [line[1:] for line in proteome.read_text().splitlines() if line.startswith('>')]
    query_headers = {fasta_accession(header.split()[0]): header for header in headers}
    if len(headers) != len(queries) or len(queries) != prior['summary']['proteome_sequences']:
        raise ValueError('Proteome duplicate identifier or inventory mismatch')
    database, hit_path = raw / 'references', raw / 'hits.tsv'
    subprocess.run(['diamond', 'makedb', '--in', str(reference_path), '--db', str(database)], check=True)
    command = ['diamond', 'blastp', '--query', str(proteome), '--db', str(database), '--out', str(hit_path),
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore', 'qcovhsp', 'scovhsp',
        '--evalue', '1e-5', '--max-target-seqs', '0', '--max-hsps', '1', '--sensitive', '--threads', '4']
    subprocess.run(command, check=True)
    hits = parse_hits(hit_path.read_text(), queries.keys(), references.keys())
    for alignments in hits.values():
        for hit in alignments:
            if hit['query_length'] != len(queries[hit['cannabis_accession']]) or hit['reference_length'] != len(references[hit['reference_accession']]['sequence']):
                raise ValueError('Alignment length differs from searched sequence')
            if any(not math.isfinite(hit[k]) or not 0 <= hit[k] <= 100 for k in ['identity_percent', 'query_coverage_percent', 'reference_coverage_percent']):
                raise ValueError('Invalid alignment percentages')
    rows, alignments = annotate(discovery, references, hits)
    candidates = sorted({p for r in rows for p in r['screened_cannabis_proteins']})
    result = {'schema': 'cannabis-carbon.phase1-new-protein-search.v1',
        'generated_at': datetime.now(timezone.utc).isoformat(),
        'source_discovery': str(source), 'source_discovery_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'proteome_path': str(proteome), 'proteome_sha256': hashlib.sha256(proteome.read_bytes()).hexdigest(),
        'reference_path': str(reference_path), 'reference_sha256': hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        'hits_path': str(hit_path), 'hits_sha256': hashlib.sha256(hit_path.read_bytes()).hexdigest(),
        'diamond_command': command, 'diamond_version': subprocess.check_output(['diamond', 'version'], text=True).strip(),
        'screen': {'minimum_identity_percent': 30, 'minimum_query_coverage_percent': 50, 'minimum_reference_coverage_percent': 50, 'maximum_evalue': 1e-5},
        'summary': {'proteome_sequences': len(queries), 'requested_references': len(ids), 'retrieved_references': len(references),
            'failed_retrieval_batches': sum(r['status'] != 'retrieved' for r in retrievals),
            'equation_gaps': len(rows), 'raw_alignments': sum(map(len, hits.values())), 'passing_alignments': len(alignments),
            'equations_with_screened_candidates': sum(bool(r['screened_cannabis_proteins']) for r in rows),
            'priority_targets_with_screened_candidates': len({t for r in rows if r['screened_cannabis_proteins'] for t in r['priority_target_ids']}),
            'distinct_cannabis_candidates': len(candidates),
            'protein_reaction_hypotheses': sum(len(r['screened_cannabis_proteins']) for r in rows)},
        'retrievals': retrievals, 'reference_sequences': list(references.values()), 'passing_alignments': alignments, 'rows': rows,
        'cannabis_candidates': [{'accession': p, 'source_header': query_headers[p], 'sequence': queries[p],
            'sequence_sha256': hashlib.sha256(queries[p].encode()).hexdigest(), 'source_url': f'https://www.uniprot.org/uniprotkb/{p}/entry'} for p in candidates],
        'claim_boundary': 'Whole-proteome homology screening, not experimental enzyme confirmation. Only passing alignments are embedded; the full alignment output is checksummed and per-equation raw counts retain weak-hit outcomes. Exact specificity, physiological direction and full CO2 pathways remain unestablished. Atom tracing is deferred.'}
    Path('data/reports/phase1-new-protein-search.json').write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(json.dumps(result['summary']), flush=True)


def export_table(report_path, output):
    raw = report_path.read_bytes()
    report = json.loads(raw)
    digest = hashlib.sha256(raw).hexdigest()
    collections = [('equation_gap', 'rows', 'reaction_id'), ('reference_sequence', 'reference_sequences', 'accession'),
        ('cannabis_candidate', 'cannabis_candidates', 'accession'), ('passing_alignment', 'passing_alignments', 'id'),
        ('retrieval', 'retrievals', 'url')]
    count = 0
    with output.open('w') as stream:
        def write(kind, identifier, row):
            stream.write(json.dumps({'record_kind': kind, 'record_id': identifier,
                'status': row.get('search_status') or row.get('status'),
                'record_json': json.dumps(row, separators=(',', ':')), 'report_sha256': digest}, separators=(',', ':')) + '\n')
        for kind, collection, key in collections:
            for row in report[collection]:
                write(kind, row[key], row); count += 1
        write('metadata', 'report', {k: v for k, v in report.items() if k not in {c for _, c, _ in collections}})
    return count + 1


if __name__ == '__main__':
    run()
