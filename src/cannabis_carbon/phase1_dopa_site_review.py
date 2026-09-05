"""Alignment-based site hypotheses, never transferred activity assignments."""
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path
from .genome import _fasta


def map_sites(fields, query, reference, positions):
    qs, qe, ss, se = map(int, fields[10:14])
    qa, sa = fields[14:16]
    if len(qa) != len(sa) or qa.replace('-', '') != query[qs-1:qe] or sa.replace('-', '') != reference[ss-1:se]:
        raise ValueError('Aligned sequences do not match pinned sequence slices')
    qpos, spos = qs - 1, ss - 1
    mapped = {}
    for q, s in zip(qa, sa):
        qpos += q != '-'
        spos += s != '-'
        if s != '-' and spos in positions:
            mapped[spos] = {'reference_position': spos, 'reference_residue': s,
                'query_position': qpos if q != '-' else None,
                'query_residue': q if q != '-' else None,
                'status': 'aligned-residue' if q != '-' else 'query-gap'}
    return [mapped.get(p, {'reference_position': p, 'reference_residue': reference[p-1],
        'query_position': None, 'query_residue': None, 'status': 'outside-local-alignment'}) for p in sorted(positions)]


def run():
    source = Path('data/reports/phase1-dopa-lyase-search.json')
    annotation_path = Path('data/raw/phase1-dopa-lyase-search/Q3IWB0.json')
    search = json.loads(source.read_text())
    annotation = json.loads(annotation_path.read_text())
    reference = annotation['sequence']['value']
    assert reference == next(r['sequence'] for r in search['reference_sequences'] if r['accession'] == 'Q3IWB0')
    features = [f for f in annotation['features'] if f['type'] in ('Binding site', 'Active site', 'Mutagenesis', 'Modified residue')]
    positions = {p for f in features for p in range(f['location']['start']['value'], f['location']['end']['value'] + 1)}
    raw = Path('data/raw/dopa-site-review')
    raw.mkdir(parents=True, exist_ok=True)
    coordinates = raw / 'coordinate-hits.tsv'
    command = list(search['diamond_command'])
    command[command.index('--out') + 1] = str(coordinates)
    i = command.index('--evalue')
    command[i:i] = ['qstart', 'qend', 'sstart', 'send', 'qseq_gapped', 'sseq_gapped']
    subprocess.run(command, check=True)
    lines = [line.split('\t') for line in coordinates.read_text().splitlines()]
    if any(len(f) != 16 for f in lines) or Counter('\t'.join(f[:10]) for f in lines) != Counter(Path(search['hits_path']).read_text().splitlines()):
        raise ValueError('Full-proteome replay differs from original screen')
    sequences = _fasta(Path(search['proteome_path']))
    leads = {a['cannabis_accession']: a['id'] for a in search['passing_alignments']}
    rows = []
    for f in lines:
        acc = f[0].split('|')[1]
        sites = map_sites(f, sequences[acc], reference, positions)
        if acc in leads:
            rows.append({'accession': acc, 'alignment_id': leads[acc], 'alignment_columns': f,
                'sites': sites, 'model_eligible': False})
    if Counter(r['accession'] for r in rows) != Counter({acc: 1 for acc in leads}):
        raise ValueError('Missing or duplicated lead')
    paths = [source, annotation_path, coordinates, Path(search['hits_path']), Path(search['proteome_path'])]
    report = {'schema': 'cannabis-dopa-site-review-v1', 'model_eligible': False,
        'reference_accession': 'Q3IWB0', 'reference_features': features,
        'coordinate_command': command, 'coordinate_replay_alignment_count': len(lines),
        'rows': sorted(rows, key=lambda r: r['accession']),
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'protein_leads': len(rows), 'reference_H89_aligned_residues': dict(Counter(next(s['query_residue'] for s in r['sites'] if s['reference_position'] == 89) or 'unmapped' for r in rows)), 'new_exact_enzyme_assignments': 0},
        'review_decision': 'Retain experimental leads only. Compare mapped site 89 with the preserved reference H89F mutagenesis annotation; that annotation concerns reference Tyr/Phe specificity, not Cannabis L-DOPA activity. Test each candidate on L-DOPA, L-tyrosine and L-phenylalanine with authentic product standards and no-enzyme, inactive-enzyme and substrate-autoxidation controls. Independently validate the reference L-DOPA activity reported without shown assay data.',
        'claim_boundary': 'Coordinates use pinned full-length UniProt sequences, not PDB construct numbering. One local sequence alignment supplies positional hypotheses, not structural equivalence or validated specificity. Reference feature evidence codes (including absent evidence) are preserved. Residue agreement or disagreement neither proves nor excludes Cannabis L-DOPA turnover. No model promotion or atom-tracing claim.'}
    Path('data/reports/phase1-dopa-site-review.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
