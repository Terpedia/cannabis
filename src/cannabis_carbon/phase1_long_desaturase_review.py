"""Review long desaturase-domain leads without overriding sequence-screen failures."""
import hashlib
import json
from pathlib import Path
from .genome import _fasta
from .phase1_family_search import parse_hits


def build():
    proteome = Path('data/raw/UP000583929.fasta')
    sequences = _fasta(proteome)
    paths = [proteome]
    searches = []
    for name in ('phase1-remaining-gap-search', 'phase1-weak-nonplant-search'):
        path = Path('data/reports/' + name + '.json')
        report = json.loads(path.read_text())
        hits_path = Path(report['hits_path'])
        if hashlib.sha256(hits_path.read_bytes()).hexdigest() != report['hits_sha256']:
            raise ValueError('Prior alignment snapshot changed')
        if hashlib.sha256(proteome.read_bytes()).hexdigest() != report['proteome_sha256']:
            raise ValueError('Proteome snapshot changed')
        hits = parse_hits(hits_path.read_text(), sequences.keys(), {r['accession'] for r in report['reference_sequences']})
        paths.extend((path, hits_path))
        searches.append((name, hits))
    rows = []
    for accession in ('A0A7J6DP00', 'A0A7J6F905'):
        path = Path('data/raw/weak-hit-domain-review/' + accession + '.json')
        protein = json.loads(path.read_text())
        if protein['sequence']['value'] != sequences[accession]:
            raise ValueError('Annotation sequence differs from pinned proteome')
        paths.append(path)
        rows.append({'accession': accession, 'sequence_length': len(sequences[accession]),
            'source_url': f'https://rest.uniprot.org/uniprotkb/{accession}.json',
            'features': protein['features'], 'comments': protein.get('comments', []),
            'domain_crossrefs': [r for r in protein['uniProtKBCrossReferences'] if r['database'] in ('Pfam', 'InterPro')],
            'prior_alignments': [{'source_report': name, 'alignment': hit}
                for name, hits in searches for group in hits.values() for hit in group
                if hit['cannabis_accession'] == accession],
            'model_eligible': False,
            'evidence_class': 'desaturase-domain-annotation-specificity-unverified-lead'})
    return {'schema': 'cannabis-long-desaturase-review-v1',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'rows': rows,
        'claim_boundary': 'Two selected longer partial-domain hits, not a new genome-wide search. Pfam FA_desaturase plus cytochrome b5 annotations motivate catalytic-domain follow-up, not exact delta-6 acyl-CoA activity. All prior alignments and their original thresholds remain unchanged; no route promotion.',
        'next_test': 'Use experimentally characterized full-length desaturases and their catalytic regions for a separate domain-aware whole-proteome search. Compare alternative substrate and double-bond positional specificities; confirm product geometry and carrier identity using authentic standards.'}


if __name__ == '__main__':
    Path('data/reports/phase1-long-desaturase-review.json').write_text(json.dumps(build(), separators=(',', ':')) + '\n')
