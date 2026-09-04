"""Proteome-wide homology search of direction-unresolved Rhea-family leads."""
import hashlib
import json
import subprocess
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from .genome import _fasta
from .phase1_search import fasta_accession


def parse_references(data, expected):
    if data.strip() and not data.lstrip().startswith(b'>'):
        raise ValueError('Expected a FASTA response')
    records = {}
    for block in data.decode('ascii').strip().split('>'):
        if not block.strip():
            continue
        lines = block.splitlines()
        accession = fasta_accession(lines[0].split()[0])
        sequence = ''.join(line.strip() for line in lines[1:])
        if accession not in expected or accession in records:
            raise ValueError('Unexpected or duplicate reference accession')
        if not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO*'):
            raise ValueError('Invalid reference sequence')
        records[accession] = {'accession': accession, 'header': lines[0], 'sequence': sequence,
                              'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest()}
    return records


def parse_hits(text, query_ids, reference_ids):
    hits = {}
    for line in text.splitlines():
        f = line.split('\t')
        if len(f) != 10:
            raise ValueError('Unexpected alignment columns')
        query, ref = fasta_accession(f[0]), fasta_accession(f[1])
        if query not in query_ids or ref not in reference_ids:
            raise ValueError('Alignment identifier absent from searched sequences')
        hit = {'cannabis_accession': query, 'reference_accession': ref,
               'identity_percent': float(f[2]), 'alignment_length': int(f[3]),
               'query_length': int(f[4]), 'reference_length': int(f[5]),
               'evalue': float(f[6]), 'bitscore': float(f[7]),
               'query_coverage_percent': float(f[8]), 'reference_coverage_percent': float(f[9])}
        hit['passes_screen'] = (hit['identity_percent'] >= 30 and hit['query_coverage_percent'] >= 50
                                and hit['reference_coverage_percent'] >= 50 and 0 <= hit['evalue'] <= 1e-5)
        hits.setdefault(ref, []).append(hit)
    return hits


def annotate_rows(discovery_rows, references, hits):
    rows = []
    for row in discovery_rows:
        refs = {r['accession'] for r in row['family_reference_annotations']}
        present = refs & references.keys()
        alignments = [hit for ref in sorted(present) for hit in hits.get(ref, [])]
        proteins = sorted({h['cannabis_accession'] for h in alignments if h['passes_screen']})
        rows.append({**row, 'sequence_search_status': 'hits-found' if alignments else 'no-hits' if present else 'no-reference-sequence',
                     'reference_sequences_present': sorted(present), 'reference_sequences_missing': sorted(refs - present),
                     'sequence_hits': alignments, 'screened_cannabis_proteins': proteins,
                     'direction_status': 'unresolved-for-requested-reaction-direction',
                     'evidence_class': 'direction-unresolved-family-homology-candidate' if proteins else 'unresolved',
                     'validation_blockers': ['requested-reaction-direction-unverified', 'exact-substrate-specificity-unverified',
                                             'catalytic-residues-and-domains-not-reviewed', 'compartment-and-expression-unverified'],
                     'proposed_test': 'Review the reference annotation and full source reaction; test candidate activity with exact substrates and cofactors in both proposed directions, including negative controls and a characterized reference control. Confirm product identity against standards before claiming the requested direction.',
                     'claim_boundary': 'Homology to a protein annotated to another direction in the same explicit Rhea family is not evidence for the requested physiological direction, exact substrate specificity, or a complete Cannabis pathway.'})
    return rows


def run():
    source = Path('data/reports/phase1-reference-discovery.json')
    discovery = json.loads(source.read_text())
    if hashlib.sha256(Path(discovery['source_overlay']).read_bytes()).hexdigest() != discovery['source_overlay_sha256']:
        raise ValueError('Discovery input overlay has changed; regenerate or resolve the snapshot before searching')
    ids = sorted({p['accession'] for r in discovery['rows'] for p in r['family_reference_annotations']})
    raw = Path('data/raw/phase1-family-search')
    raw.mkdir(parents=True, exist_ok=True)
    references, retrievals = {}, []
    for start in range(0, len(ids), 60):
        batch = ids[start:start + 60]
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({
            'query': ' OR '.join('accession:' + accession for accession in batch), 'format': 'fasta'})
        retrieval = {'requested_accessions': batch, 'url': url}
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                data = response.read()
            parsed = parse_references(data, set(batch))
            digest = hashlib.sha256(data).hexdigest()
            snapshot = raw / (digest + '.fasta')
            snapshot.write_bytes(data)
            retrieval.update(status='retrieved', sha256=digest, snapshot=str(snapshot), missing_accessions=sorted(set(batch) - parsed.keys()))
            references.update(parsed)
        except (OSError, ValueError) as error:
            retrieval.update(status='retrieval-failed', reason=str(error))
        retrievals.append(retrieval)
        print(f'References {min(start + 60, len(ids))}/{len(ids)}: {retrieval["status"]}', flush=True)
    if not references:
        raise ValueError('No reference sequences retrieved; no search performed')
    reference_path = raw / 'references.fasta'
    reference_path.write_text(''.join(f'>{accession} {record["header"]}\n{record["sequence"]}\n' for accession, record in sorted(references.items())))
    proteome = Path('data/raw/UP000583929.fasta')
    queries = _fasta(proteome)
    query_headers = {fasta_accession(line[1:].split()[0]): line[1:] for line in proteome.read_text().splitlines() if line.startswith('>')}
    database, output = raw / 'references', raw / 'hits.tsv'
    subprocess.run(['diamond', 'makedb', '--in', str(reference_path), '--db', str(database)], check=True)
    command = ['diamond', 'blastp', '--query', str(proteome), '--db', str(database), '--out', str(output),
               '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore', 'qcovhsp', 'scovhsp',
               '--evalue', '1e-5', '--max-target-seqs', '0', '--max-hsps', '1', '--sensitive', '--threads', '4']
    subprocess.run(command, check=True)
    hits = parse_hits(output.read_text(), queries.keys(), references.keys())
    rows = annotate_rows(discovery['rows'], references, hits)
    candidate_ids = sorted({p for r in rows for p in r['screened_cannabis_proteins']})
    result = {'schema': 'cannabis-carbon.phase1-family-protein-search.v1',
              'generated_at': datetime.now(timezone.utc).isoformat(),
              'source_discovery': str(source), 'source_discovery_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'proteome_path': str(proteome), 'proteome_sha256': hashlib.sha256(proteome.read_bytes()).hexdigest(),
              'reference_sha256': hashlib.sha256(reference_path.read_bytes()).hexdigest(),
              'hits_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
              'diamond_command': command, 'diamond_version': subprocess.check_output(['diamond', 'version'], text=True).strip(),
              'retrievals': retrievals, 'reference_sequences': list(references.values()),
              'cannabis_candidates': [{'accession': p, 'source_header': query_headers[p], 'sequence': queries[p],
                                       'sequence_sha256': hashlib.sha256(queries[p].encode()).hexdigest(),
                                       'source_url': f'https://www.uniprot.org/uniprotkb/{p}/entry'} for p in candidate_ids],
              'summary': {'proteome_sequences': len(queries), 'requested_references': len(ids), 'retrieved_references': len(references),
                          'gap_variants': len(rows), 'variants_with_hits': sum(bool(r['sequence_hits']) for r in rows),
                          'variants_with_screened_family_candidates': sum(bool(r['screened_cannabis_proteins']) for r in rows),
                          'distinct_screened_cannabis_proteins': len({p for r in rows for p in r['screened_cannabis_proteins']}),
                          'variants_with_weak_hits_only': sum(bool(r['sequence_hits']) and not r['screened_cannabis_proteins'] for r in rows),
                          'variants_with_no_hits': sum(r['sequence_search_status'] == 'no-hits' for r in rows),
                          'variants_without_reference_sequences': sum(r['sequence_search_status'] == 'no-reference-sequence' for r in rows),
                          'protein_reaction_hypotheses': sum(len(r['screened_cannabis_proteins']) for r in rows)},
              'rows': rows, 'claim_boundary': 'Direction-unresolved family homology is a separate evidence layer and must not silently increase direction-specific or confirmed-enzyme completeness.'}
    Path('data/reports/phase1-family-protein-search.json').write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(result['summary'])


def export_table(report_path, output):
    data = report_path.read_bytes()
    report = json.loads(data)
    rows = []
    for row in report['rows']:
        rows.append({'reaction_id': row['reaction_id'], 'reaction_smarts': row['reaction_smarts'],
                     'sequence_search_status': row['sequence_search_status'], 'direction_status': row['direction_status'],
                     'screened_protein_count': len(row['screened_cannabis_proteins']),
                     'candidate_proteins_json': json.dumps(row['screened_cannabis_proteins']),
                     'reference_annotations_json': json.dumps(row['family_reference_annotations']),
                     'sequence_hits_json': json.dumps(row['sequence_hits']),
                     'reference_sequences_missing_json': json.dumps(row['reference_sequences_missing']),
                     'validation_blockers_json': json.dumps(row['validation_blockers']),
                     'proposed_test': row['proposed_test'], 'claim_boundary': row['claim_boundary'],
                     'report_sha256': hashlib.sha256(data).hexdigest(),
                     'source_discovery_sha256': report['source_discovery_sha256'],
                     'proteome_sha256': report['proteome_sha256'], 'hits_sha256': report['hits_sha256']})
    output.write_text(''.join(json.dumps(row, separators=(',', ':')) + '\n' for row in rows))
    return len(rows)


if __name__ == '__main__':
    run()
