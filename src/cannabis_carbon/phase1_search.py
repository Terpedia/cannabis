"""Search the Cannabis proteome against reaction-linked UniProt references."""
import hashlib
import json
import subprocess
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from .genome import _fasta


def run():
    queue_path = Path('data/reports/phase1-enzyme-discovery-queue.json')
    queue = json.loads(queue_path.read_text())['rows']
    ids = sorted({i for r in queue if r['balance_status'] == 'balanced' for i in r['source_uniprot_ids']})
    # UniParc identifiers need a separate service; retain them as missing
    # references instead of sending an invalid UniProtKB accession query.
    kb_ids = [i for i in ids if not i.startswith('UPI')]
    query = ' OR '.join('accession:' + i for i in kb_ids)
    url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query': query, 'format': 'fasta'})
    reference = Path('data/raw/phase1-reaction-references.fasta')
    with urllib.request.urlopen(url, timeout=60) as response:
        reference.write_bytes(response.read())
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
        ref, cannabis = f[1].split('|')[1], f[0].split('|')[1]
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
