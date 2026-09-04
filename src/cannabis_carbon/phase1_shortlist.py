"""Reproducible protein/reaction shortlist for Phase 1 experimental review."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from .genome import _fasta
from .phase1_search import fasta_accession


def build_shortlist(search_path, proteome_path, output):
    search = json.loads(search_path.read_text())
    digest = hashlib.sha256(proteome_path.read_bytes()).hexdigest()
    if digest != search['proteome_sha256']:
        raise ValueError('Proteome differs from the searched sequence snapshot')
    sequences = _fasta(proteome_path)
    headers = {fasta_accession(line[1:].split()[0]): line[1:]
               for line in proteome_path.read_text().splitlines() if line.startswith('>')}
    proteins = defaultdict(list)
    for reaction in search['rows']:
        grouped = defaultdict(list)
        for hit in reaction['sequence_hits']:
            if hit['passes_screen']:
                grouped[hit['cannabis_accession']].append(hit)
        for accession, hits in grouped.items():
            best = min(hits, key=lambda h: (h['evalue'], -h['bitscore'], h['reference_accession']))
            proteins[accession].append({
                'reaction_id': reaction['reaction_id'], 'reaction_smarts': reaction['reaction_smarts'],
                'balance_status': reaction['balance_status'], 'source_urls': reaction['source_urls'],
                'candidate_cannabisdb_ids': reaction['candidate_cannabisdb_ids'],
                'best_alignment': best, 'reference_alignments': hits,
                'status': 'homology_candidate',
                'validation_blockers': ['reaction-specific-catalysis-unverified', 'physiological-direction-unverified',
                                        'catalytic-residues-and-domains-not-reviewed', 'compartment-and-expression-unverified'],
                'proposed_test': 'Review source reaction direction and full stoichiometry; assay the candidate with the exact substrates and required cofactors, appropriate negative controls and a characterized positive control; confirm product identity against standards. Test competing substrate/product hypotheses for the same protein.',
            })
    rows = []
    for accession, hypotheses in sorted(proteins.items()):
        hypotheses.sort(key=lambda h: (h['best_alignment']['evalue'], -h['best_alignment']['bitscore'], h['reaction_id'], h['reaction_smarts']))
        sequence = sequences[accession]
        rows.append({'accession': accession, 'source_header': headers[accession],
                     'source_url': f'https://www.uniprot.org/uniprotkb/{accession}/entry',
                     'sequence': sequence, 'sequence_length': len(sequence),
                     'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                     'reaction_hypothesis_count': len(hypotheses), 'hypotheses': hypotheses})
    result = {'schema': 'cannabis-carbon.phase1-experimental-shortlist.v1',
              'source_search': str(search_path), 'search_sha256': hashlib.sha256(search_path.read_bytes()).hexdigest(),
              'proteome_sha256': digest, 'protein_count': len(rows),
              'protein_reaction_hypothesis_count': sum(r['reaction_hypothesis_count'] for r in rows),
              'ranking': 'Within each protein: E-value ascending, bitscore descending; ranking is alignment evidence, not probability of catalytic activity.',
              'proteins': rows,
              'claim_boundary': 'FASTA headers retain source annotations, which may describe characterized activities. Proposed reaction variants remain unverified; related variants and repeated reference hits are not independent confirmations.'}
    output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    return {k: result[k] for k in ('protein_count', 'protein_reaction_hypothesis_count')}


if __name__ == '__main__':
    print(build_shortlist(Path('data/reports/phase1-targeted-protein-search.json'),
                         Path('data/raw/UP000583929.fasta'), Path('data/reports/phase1-experimental-shortlist.json')))
