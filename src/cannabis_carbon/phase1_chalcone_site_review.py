"""Replay complete CHI screen and map the characterized reference's sites."""
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path
from .genome import _fasta
from .phase1_dopa_site_review import map_sites

REFERENCE = 'P28012'


def run():
    source = Path('data/reports/phase1-chalcone-search.json')
    annotation_path = Path('data/raw/chalcone-annotations/P28012.json')
    search = json.loads(source.read_text())
    annotation = json.loads(annotation_path.read_text())
    references = {r['accession']: r['sequence'] for r in search['reference_sequences']}
    if annotation['sequence']['value'] != references[REFERENCE]:
        raise ValueError('Reference annotation sequence mismatch')
    features = [f for f in annotation['features'] if f['type'] in ('Binding site', 'Active site', 'Site', 'Mutagenesis')]
    positions = {p for f in features for p in range(f['location']['start']['value'], f['location']['end']['value'] + 1)}
    raw = Path('data/raw/chalcone-site-review')
    raw.mkdir(parents=True, exist_ok=True)
    coordinates = raw / 'coordinate-hits.tsv'
    command = list(search['diamond_command'])
    command[command.index('--out') + 1] = str(coordinates)
    i = command.index('--evalue')
    command[i:i] = ['qstart', 'qend', 'sstart', 'send', 'qseq_gapped', 'sseq_gapped']
    subprocess.run(command, check=True)
    lines = [line.split('\t') for line in coordinates.read_text().splitlines()]
    if any(len(f) != 16 for f in lines) or Counter('\t'.join(f[:10]) for f in lines) != Counter(Path(search['hits_path']).read_text().splitlines()):
        raise ValueError('Full-proteome coordinate replay changed original alignments')
    queries = _fasta(Path(search['proteome_path']))
    passing = {(a['cannabis_accession'], a['reference_accession']): a['id'] for a in search['passing_alignments']}
    mapped = []
    for f in lines:
        acc, ref = f[0].split('|')[1], f[1]
        sites = map_sites(f, queries[acc], references[ref], positions if ref == REFERENCE else set())
        if ref == REFERENCE:
            mapped.append({'accession': acc, 'alignment_columns': f, 'sites': sites,
                'passing_alignment_id': passing.get((acc, ref)), 'passes_original_screen': (acc, ref) in passing,
                'model_eligible': False})
    absent = sorted({p['accession'] for p in search['cannabis_candidates']} - {r['accession'] for r in mapped})
    paths = [source, annotation_path, coordinates, Path(search['hits_path']), Path(search['proteome_path'])]
    report = {'schema': 'cannabis-chalcone-site-review-v1', 'reference_accession': REFERENCE,
        'reference_features': features, 'model_eligible': False, 'coordinate_command': command,
        'coordinate_replay_alignment_count': len(lines), 'rows': mapped,
        'other_reference_leads_without_P28012_alignment': absent,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'reference_site_positions': sorted(positions), 'P28012_alignments': len(mapped),
            'P28012_passing_alignments': sum(r['passes_original_screen'] for r in mapped), 'new_exact_enzyme_assignments': 0},
        'claim_boundary': 'All original hits are replayed, including failed thresholds. Sites use pinned UniProt numbering, not assumed PDB construct numbering. No reported local alignment is recorded explicitly, not interpreted as enzyme absence. Sequence alignment suggests positional correspondence, not structural equivalence, catalysis or stereoselectivity. Residue conservation cannot by itself promote a reference-class lead to exact Cannabis reaction evidence. Atom tracing and candidate model unchanged.'}
    Path('data/reports/phase1-chalcone-site-review.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))
    print(json.dumps([{'accession': r['accession'], 'passes': r['passes_original_screen'], 'sites': r['sites']} for r in mapped]))


if __name__ == '__main__':
    run()
