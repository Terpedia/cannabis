"""Replay a pinned predicted CDS from genomic bases, retaining every codon."""
import hashlib
import itertools
import json
import urllib.parse
import urllib.request
from pathlib import Path

# NCBI standard nuclear genetic code, codons ordered T,C,A,G at each position.
CODE = dict(zip((''.join(c) for c in itertools.product('TCAG', repeat=3)),
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


def reconstruct(sequence, region_start, exons):
    bases, coords = [], []
    for exon in exons:
        if exon['strand'] != -1:
            raise ValueError('Expected reverse-strand model')
        start, end = exon['start'], exon['end']
        if not region_start <= start <= end < region_start + len(sequence):
            raise ValueError('Exon outside retrieved interval')
        fragment = sequence[start-region_start:end-region_start+1]
        if set(fragment) - set('ACGT'):
            raise ValueError('Ambiguous genomic base')
        bases.append(fragment.translate(str.maketrans('ACGT', 'TGCA'))[::-1])
        coords.extend(range(end, start-1, -1))
    cds = ''.join(bases)
    if len(cds) % 3:
        raise ValueError('Incomplete codon')
    codons = [{'codon_index': i//3+1, 'codon': cds[i:i+3], 'amino_acid': CODE[cds[i:i+3]],
        'genomic_positions': coords[i:i+3]} for i in range(0, len(cds), 3)]
    return cds, codons


def run():
    parent = Path('data/reports/phase1-chalcone-gene-model.json')
    model = json.loads(parent.read_text())['source_record']
    exons = model['transcript_order_exons']
    accessions = {e['accession'] for e in exons}
    if len(accessions) != 1:
        raise ValueError('Multiple genomic accessions')
    accession = next(iter(accessions))
    start, end = min(e['start'] for e in exons), max(e['end'] for e in exons)
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urllib.parse.urlencode({
        'db': 'nuccore', 'id': accession, 'rettype': 'fasta', 'retmode': 'text',
        'seq_start': start, 'seq_stop': end, 'strand': 1})
    path = Path('data/raw/chalcone-genomic-translation') / f'{accession}-{start}-{end}.fasta'
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = response.read()
        if not payload.startswith(('>' + accession + ':' + str(start) + '-' + str(end) + ' ').encode()):
            raise ValueError('Unexpected genomic response identity')
        path.write_bytes(payload)
    lines = path.read_text().splitlines()
    if not lines[0].startswith(f'>{accession}:{start}-{end} '):
        raise ValueError('Cached genomic identity mismatch')
    sequence = ''.join(lines[1:]).upper()
    if len(sequence) != end-start+1:
        raise ValueError('Truncated genomic sequence')
    cds, codons = reconstruct(sequence, start, exons)
    translation = ''.join(c['amino_acid'] for c in codons)
    if translation != model['sequence'] + '*' or not cds.startswith('ATG'):
        raise ValueError('Predicted protein and spliced genomic translation differ')
    report = {'schema': 'cannabis-chalcone-genomic-translation-v1', 'model_eligible': False,
        'genetic_code': 1, 'source_url': url, 'genomic_accession': accession,
        'genomic_region': {'start': start, 'end': end, 'retrieved_strand': 1},
        'transcript_order_exons': exons, 'spliced_cds': cds, 'codons': codons,
        'translation_including_stop': translation,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (parent, path)},
        'summary': {'genomic_bases': len(sequence), 'coding_bases': len(cds), 'protein_residues': len(translation)-1,
            'terminal_stop': codons[-1]['codon'], 'exact_predicted_protein_match': True},
        'claim_boundary': 'Exact genomic-to-predicted-protein replay, including codon-level genomic coordinates and terminal stop. This validates sequence arithmetic for the supplied six-exon model, not its biological correctness, splicing, expression, mature termini, or enzymatic function. No independent transcript evidence, model promotion, or carbon-atom mapping.'}
    Path('data/reports/phase1-chalcone-genomic-translation.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
