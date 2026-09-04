"""Search the Cannabis proteome against reaction-linked UniProt/NCBI references."""
import hashlib
import json
import subprocess
import urllib.parse
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

from .genome import _fasta


def fasta_accession(identifier):
    parts = identifier.split('|')
    return parts[1] if len(parts) >= 3 else identifier.split()[0]


def validate_single_reference(accession, data):
    """Require the exact versioned accession and a single nonempty sequence."""
    lines = data.decode('ascii').strip().splitlines()
    if not lines or not lines[0].startswith('>' + accession + ' '):
        raise ValueError('Returned FASTA identifier does not match request')
    if sum(line.startswith('>') for line in lines) != 1:
        raise ValueError('Expected exactly one reference sequence')
    sequence = ''.join(line.strip() for line in lines[1:])
    if not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO*'):
        raise ValueError('Invalid or empty protein sequence')
    return lines[0][1:]


def annotated_cds_reference(accession, data):
    """Use a single source-annotated CDS, never invent a nucleotide translation."""
    records = ET.fromstring(data).findall('GBSeq')
    if len(records) != 1:
        raise ValueError('Expected exactly one nucleotide record')
    record = records[0]
    version = record.findtext('GBSeq_accession-version')
    expected = version if '.' in accession else record.findtext('GBSeq_primary-accession')
    if accession != expected:
        raise ValueError('Nucleotide accession does not match request')
    cds = [f for f in record.findall('./GBSeq_feature-table/GBFeature') if f.findtext('GBFeature_key') == 'CDS']
    if len(cds) != 1:
        raise ValueError('Source CDS identity is ambiguous or missing')
    qualifiers = {q.findtext('GBQualifier_name'): q.findtext('GBQualifier_value')
                  for q in cds[0].findall('./GBFeature_quals/GBQualifier')}
    protein, sequence = qualifiers.get('protein_id'), qualifiers.get('translation')
    if not protein or not sequence:
        raise ValueError('CDS lacks a protein accession or annotated translation')
    organism = record.findtext('GBSeq_organism')
    header = f"{protein} {qualifiers.get('product', 'annotated CDS')} [{organism}] nucleotide={version}"
    fasta = f'>{header}\n{sequence}\n'.encode('ascii')
    validate_single_reference(protein, fasta)
    return fasta, {'protein_accession': protein, 'nucleotide_accession_version': version,
                   'source_header': header, 'cds_location': cds[0].findtext('GBFeature_location'),
                   'pmids': [p.text for p in record.findall('./GBSeq_references/GBReference/GBReference_pubmed')],
                   'translation_method': 'NCBI-source-annotated-single-CDS; not independently translated'}


def run():
    queue_path = Path('data/reports/phase1-enzyme-discovery-queue.json')
    queue = json.loads(queue_path.read_text())['rows']
    ids = sorted({i for r in queue if r['balance_status'] == 'balanced' for i in r['source_uniprot_ids']})
    genbank_ids = sorted({i for r in queue if r['balance_status'] == 'balanced' for i in r['source_genbank_ids']})
    # UniParc records retain their archive ID as the alignment reference ID.
    kb_ids = [i for i in ids if not i.startswith('UPI')]
    query = ' OR '.join('accession:' + i for i in kb_ids)
    url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query': query, 'format': 'fasta'})
    reference = Path('data/raw/phase1-reaction-references.fasta')
    with urllib.request.urlopen(url, timeout=60) as response:
        reference_data = response.read()
    retrievals = []
    for accession in (i for i in ids if i.startswith('UPI')):
        archive_url = f'https://rest.uniprot.org/uniparc/{accession}.fasta'
        try:
            with urllib.request.urlopen(archive_url, timeout=30) as response:
                fasta = response.read()
            if not fasta.startswith(('>' + accession + ' ').encode()):
                raise ValueError('Returned FASTA identifier does not match request')
            reference_data += b'\n' + fasta
            retrievals.append({'accession': accession, 'url': archive_url, 'status': 'retrieved',
                               'sha256': hashlib.sha256(fasta).hexdigest()})
        except (urllib.error.URLError, TimeoutError, ValueError) as error:
            retrievals.append({'accession': accession, 'url': archive_url, 'status': 'unresolved', 'reason': str(error)})
    genbank_retrievals = []
    resolved_ids = {accession: accession for accession in ids}
    for accession in genbank_ids:
        ncbi_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urllib.parse.urlencode(
            {'db': 'protein', 'id': accession, 'rettype': 'fasta', 'retmode': 'text'})
        try:
            with urllib.request.urlopen(ncbi_url, timeout=30) as response:
                fasta = response.read()
            header = validate_single_reference(accession, fasta)
            reference_data += b'\n' + fasta
            resolved_ids[accession] = accession
            genbank_retrievals.append({'accession': accession, 'url': ncbi_url, 'status': 'retrieved',
                                       'source_header': header, 'sha256': hashlib.sha256(fasta).hexdigest()})
        except (urllib.error.URLError, TimeoutError, ValueError) as error:
            nucleotide_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urllib.parse.urlencode(
                {'db': 'nuccore', 'id': accession, 'rettype': 'gb', 'retmode': 'xml'})
            try:
                with urllib.request.urlopen(nucleotide_url, timeout=30) as response:
                    source_xml = response.read()
                fasta, metadata = annotated_cds_reference(accession, source_xml)
                source_digest = hashlib.sha256(source_xml).hexdigest()
                snapshot = Path('data/raw/phase1-nucleotide-references') / (source_digest + '.xml')
                snapshot.parent.mkdir(parents=True, exist_ok=True)
                snapshot.write_bytes(source_xml)
                reference_data += b'\n' + fasta
                resolved_ids[accession] = metadata['protein_accession']
                genbank_retrievals.append({'accession': accession, 'url': nucleotide_url, 'status': 'retrieved-cds',
                    'protein_request_url': ncbi_url, 'protein_request_error': str(error),
                    'source_xml_sha256': source_digest, 'source_xml_snapshot': str(snapshot),
                    'sha256': hashlib.sha256(fasta).hexdigest(), **metadata})
            except (urllib.error.URLError, TimeoutError, ValueError, ET.ParseError) as cds_error:
                genbank_retrievals.append({'accession': accession, 'url': ncbi_url, 'status': 'unresolved',
                    'reason': str(error), 'nucleotide_url': nucleotide_url, 'nucleotide_reason': str(cds_error)})
    reference.write_bytes(reference_data)
    sequences = _fasta(reference)
    if not sequences:
        raise ValueError('No reference sequences returned')
    proteome = Path('data/raw/UP000583929.fasta')
    database = 'data/raw/phase1-reaction-references'
    hits_path = Path('data/raw/phase1-cannabis-hits.tsv')
    subprocess.run(['diamond', 'makedb', '--in', str(reference), '--db', database], check=True)
    command = ['diamond', 'blastp', '--query', str(proteome), '--db', database,
               '--out', str(hits_path), '--outfmt', '6', 'qseqid', 'sseqid', 'pident',
               'length', 'qlen', 'slen', 'evalue', 'bitscore', 'qcovhsp', 'scovhsp',
               '--evalue', '1e-5', '--max-target-seqs', '0', '--threads', '4']
    subprocess.run(command, check=True)
    hits = {}
    for line in hits_path.read_text().splitlines():
        f = line.split('\t')
        ref, cannabis = fasta_accession(f[1]), fasta_accession(f[0])
        hit = {'cannabis_accession': cannabis, 'reference_accession': ref,
               'identity_percent': float(f[2]), 'query_coverage_percent': float(f[8]),
               'reference_coverage_percent': float(f[9]), 'evalue': float(f[6]), 'bitscore': float(f[7])}
        hit['passes_screen'] = float(f[2]) >= 30 and float(f[8]) >= 50 and float(f[9]) >= 50
        hits.setdefault(ref, []).append(hit)
    rows = []
    for r in queue:
        if r['balance_status'] != 'balanced':
            continue
        source_refs = sorted(set(r['source_uniprot_ids']) | set(r['source_genbank_ids']))
        refs = sorted({resolved_ids[ref] for ref in source_refs if ref in resolved_ids})
        candidates = [h for ref in refs for h in hits.get(ref, [])]
        rows.append({**r, 'reference_sequences_present': sorted(set(refs) & sequences.keys()),
                     'reference_sequences_missing': [ref for ref in source_refs if resolved_ids.get(ref) not in sequences],
                     'source_reference_resolution': {ref: resolved_ids.get(ref) for ref in source_refs},
                     'sequence_hits': candidates,
                     'screened_candidate_count': len({h['cannabis_accession'] for h in candidates if h['passes_screen']}),
                     'search_status': 'hits-found' if candidates else 'no-hits' if set(refs) & sequences.keys() else 'no-reference-sequence'})
    result = {'schema': 'cannabis-carbon.phase1-targeted-protein-search.v1',
              'generated_at': datetime.now(timezone.utc).isoformat(), 'source_queue': str(queue_path),
              'reference_url': url, 'requested_reference_count': len(set(ids) | set(genbank_ids)), 'retrieved_reference_count': len(sequences),
              'requested_uniprot_reference_count': len(ids), 'requested_genbank_reference_count': len(genbank_ids),
              'missing_reference_ids': sorted(ref for ref in set(ids) | set(genbank_ids) if resolved_ids.get(ref) not in sequences),
              'uniparc_retrievals': retrievals,
              'genbank_retrievals': genbank_retrievals,
              'distinct_screened_cannabis_proteins': len({h['cannabis_accession'] for r in rows for h in r['sequence_hits'] if h['passes_screen']}),
              'proteome_sha256': hashlib.sha256(proteome.read_bytes()).hexdigest(),
              'reference_sha256': hashlib.sha256(reference.read_bytes()).hexdigest(),
              'command': command, 'diamond_version': subprocess.check_output(['diamond', 'version'], text=True).strip(),
              'balanced_reaction_variants': len(rows),
              'variants_with_hits': sum(bool(r['sequence_hits']) for r in rows),
              'variants_with_screened_candidates': sum(r['screened_candidate_count'] > 0 for r in rows),
              'rows': rows,
              'claim_boundary': 'Homology candidates require catalytic-residue, substrate-specificity, and species-specific experimental validation. Missing sequences and no-hit searches remain explicit.'}
    Path('data/reports/phase1-targeted-protein-search.json').write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print({k: v for k, v in result.items() if k.endswith('_count') or k.startswith('variants_')})


if __name__ == '__main__':
    run()
