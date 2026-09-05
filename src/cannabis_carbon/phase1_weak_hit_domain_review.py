"""Evidence-preserving review of a short accessory-domain homology hit."""
import hashlib
import json
from collections import Counter
from pathlib import Path


def build():
    source = Path('data/reports/phase1-weak-nonplant-search.json')
    search = json.loads(source.read_text())
    raw = Path('data/raw/weak-hit-domain-review')
    coordinates = raw / 'hits-with-coordinates.tsv'
    lines = [line.split('\t') for line in coordinates.read_text().splitlines()]
    original = Path(search['hits_path'])
    if Counter('\t'.join(f[:10]) for f in lines) != Counter(original.read_text().splitlines()):
        raise ValueError('Coordinate rerun differs from original alignments')
    hit = next(h for h in search['passing_alignments'] if h['cannabis_accession'] == 'A0A7J6FB06')
    proteins = []
    sources = [source, coordinates, original]
    for accession, collection in [('A0A7J6FB06', 'cannabis_candidates'), ('A0ACB9UWW7', 'reference_sequences')]:
        snapshot = raw / (accession + '.json')
        sources.append(snapshot)
        protein = json.loads(snapshot.read_text())
        expected = next(p for p in search[collection] if p['accession'] == accession)
        if protein['sequence']['value'] != expected['sequence']:
            raise ValueError('Domain snapshot sequence differs from searched sequence')
        proteins.append({'accession': accession, 'url': f'https://rest.uniprot.org/uniprotkb/{accession}.json',
                         'sequence_length': protein['sequence']['length'],
                         'features': protein.get('features', []), 'comments': protein.get('comments', []),
                         'domain_crossrefs': [r for r in protein['uniProtKBCrossReferences'] if r['database'] in ('Pfam', 'InterPro')]})
    alignment = next(f for f in lines if 'A0A7J6FB06' in f[0] and f[1] == 'A0ACB9UWW7')
    return {'schema': 'cannabis-weak-hit-domain-review-v1',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'alignment': hit, 'coordinates': dict(zip(('qstart', 'qend', 'sstart', 'send'), map(int, alignment[10:]))),
        'proteins': proteins,
        'decision': {'model_eligible': False, 'evidence_class': 'accessory-domain-match-not-catalytic-enzyme-assignment',
            'reason': 'Alignment overlaps cytochrome b5 heme-binding annotations in both sequences. The 106-residue reference carries automated ARBA desaturase annotations and a conserved-residue caution. These data do not establish a complete desaturase catalytic domain or activity in the 133-residue Cannabis protein.',
            'boundary': 'Annotation-based domain review, not biochemical characterization or proof that the reference annotation is false. The original passing screen is preserved unchanged. An electron-transfer partner role is a hypothesis, not an established interaction.',
            'next_test': 'Review full-length desaturase references and catalytic-domain alignments independently of b5 similarity. Inspect longer Cannabis partial matches A0A7J6DP00 and A0A7J6F905 without lowering whole-protein thresholds globally; assay exact substrate conversion with and without the putative b5 partner.'}}


if __name__ == '__main__':
    Path('data/reports/phase1-weak-hit-domain-review.json').write_text(json.dumps(build(), separators=(',', ':')) + '\n')
