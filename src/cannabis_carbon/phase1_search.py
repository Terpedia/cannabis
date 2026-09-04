"""Search the Cannabis proteome against reaction-linked UniProt references."""
import hashlib
import json
import subprocess
import urllib.parse
import urllib.request
import urllib.error
from datetime import datetime, timezone
from pathlib import Path

from .genome import _fasta


def fasta_accession(identifier):
    parts = identifier.split('|')
    return parts[1] if len(parts) >= 3 else identifier.split()[0]


def run():
    queue_path = Path('data/reports/phase1-enzyme-discovery-queue.json')
    queue = json.loads(queue_path.read_text())['rows']
    ids = sorted({i for r in queue if r['balance_status'] == 'balanced' for i in r['source_uniprot_ids']})
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
        refs = r['source_uniprot_ids']
        candidates = [h for ref in refs for h in hits.get(ref, [])]
        rows.append({**r, 'reference_sequences_present': sorted(set(refs) & sequences.keys()),
                     'reference_sequences_missing': sorted(set(refs) - sequences.keys()),
                     'sequence_hits': candidates,
                     'screened_candidate_count': len({h['cannabis_accession'] for h in candidates if h['passes_screen']}),
                     'search_status': 'hits-found' if candidates else 'no-hits' if set(refs) & sequences.keys() else 'no-reference-sequence'})
    result = {'schema': 'cannabis-carbon.phase1-targeted-protein-search.v1',
              'generated_at': datetime.now(timezone.utc).isoformat(), 'source_queue': str(queue_path),
              'reference_url': url, 'requested_reference_count': len(ids), 'retrieved_reference_count': len(sequences),
              'missing_reference_ids': sorted(set(ids) - sequences.keys()),
              'uniparc_retrievals': retrievals,
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
