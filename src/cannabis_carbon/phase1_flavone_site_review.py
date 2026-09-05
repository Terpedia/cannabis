"""Replay comparator alignments and preserve seven-site positional hypotheses."""
import hashlib
import json
import re
import subprocess
from collections import Counter
from pathlib import Path
from .genome import _fasta
from .phase1_dopa_site_review import map_sites


def verify_numbering(review, fht, fnsi):
    sites = {}
    for substitution in review['tested_substitutions']:
        match = re.fullmatch(r'([A-Z])(\d+)([A-Z])', substitution)
        if not match:
            raise ValueError('Invalid substitution')
        before, pos, after = match.groups()
        pos = int(pos)
        if pos in sites or not 1 <= pos <= min(len(fht), len(fnsi)):
            raise ValueError('Duplicate or out-of-range site')
        if fht[pos-1] != before or fnsi[pos-1] != after:
            raise ValueError('Reference numbering mismatch')
        sites[pos] = {'fht_residue': before, 'fnsi_residue': after}
    return sites


def restore_explicit_masks(fields, query, reference):
    """Permit reported X masking only, retaining each replacement as evidence."""
    restored, masks = list(fields), []
    for side, sequence, start_col, end_col, aligned_col in (
        ('query', query, 10, 11, 14), ('reference', reference, 12, 13, 15)):
        start, end = int(fields[start_col]), int(fields[end_col])
        if not 1 <= start <= end <= len(sequence):
            raise ValueError('Invalid alignment bounds')
        pos, chars = start - 1, []
        for char in fields[aligned_col]:
            if char == '-':
                chars.append(char)
                continue
            if pos >= end:
                raise ValueError('Alignment exceeds coordinate span')
            expected = sequence[pos]
            pos += 1
            if char != expected:
                if char != 'X':
                    raise ValueError('Non-mask sequence mismatch')
                masks.append({'side': side, 'position': pos, 'reported_residue': char, 'pinned_residue': expected})
            chars.append(expected)
        if pos != end:
            raise ValueError('Alignment shorter than coordinate span')
        restored[aligned_col] = ''.join(chars)
    return restored, masks


def run():
    source = Path('data/reports/phase1-flavone-fht-comparison.json')
    review_path = Path('data/curation/flavone-seven-site-review.json')
    fht_path = Path('data/raw/flavone-fht-comparison/Q7XZQ7.json')
    fnsi_path = Path('data/raw/flavone-annotations/Q7XZQ8.json')
    original_path = Path('data/reports/phase1-flavone-search.json')
    comparison = json.loads(source.read_text())
    review = json.loads(review_path.read_text())
    original = json.loads(original_path.read_text())
    for filename, digest in comparison['source_sha256'].items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest() != digest:
            raise ValueError('Changed comparison source')
    fht = json.loads(fht_path.read_text())['sequence']['value']
    fnsi = json.loads(fnsi_path.read_text())['sequence']['value']
    if fnsi != next(r['sequence'] for r in original['reference_sequences'] if r['accession'] == 'Q7XZQ8'):
        raise ValueError('FNS-I reference differs from original search')
    positions = verify_numbering(review, fht, fnsi)
    raw = Path('data/raw/flavone-site-review')
    raw.mkdir(parents=True, exist_ok=True)
    coordinates = raw / 'coordinate-hits.tsv'
    command = list(comparison['diamond_command'])
    prior_hits = Path(command[command.index('--out') + 1])
    command[command.index('--out') + 1] = str(coordinates)
    i = command.index('--evalue')
    command[i:i] = ['qstart', 'qend', 'sstart', 'send', 'qseq_gapped', 'sseq_gapped']
    subprocess.run(command, check=True)
    lines = [line.split('\t') for line in coordinates.read_text().splitlines()]
    if any(len(f) != 16 for f in lines) or Counter('\t'.join(f[:10]) for f in lines) != Counter(prior_hits.read_text().splitlines()):
        raise ValueError('Full alignment replay differs')
    proteome = Path(original['proteome_path'])
    queries = _fasta(proteome)
    refs = {r['accession']: r['annotation']['sequence']['value'] for r in comparison['references']}
    mapped, masked_alignments = {}, []
    for fields in lines:
        acc, ref = fields[0].split('|')[1], fields[1]
        restored, masks = restore_explicit_masks(fields, queries[acc], refs[ref])
        sites = map_sites(restored, queries[acc], refs[ref], positions if ref == 'Q7XZQ7' else [])
        if masks:
            masked_alignments.append({'accession': acc, 'reference_accession': ref, 'masks': masks})
        for site in sites:
            site['reported_masking'] = [m for m in masks if
                (m['side'] == 'query' and m['position'] == site['query_position']) or
                (m['side'] == 'reference' and m['position'] == site['reference_position'])]
        if ref == 'Q7XZQ7':
            if acc in mapped:
                raise ValueError('Duplicate FHT alignment')
            mapped[acc] = {'alignment_columns': fields, 'reported_masks': masks,
                'sites': [dict(s, **positions[s['reference_position']]) for s in sites]}
    rows = []
    for lead in comparison['rows']:
        acc = lead['accession']
        hit = next(c for c in lead['comparators'] if c['reference_accession'] == 'Q7XZQ7')
        rows.append({'accession': acc, 'model_eligible': False, 'comparison_status': hit['status'],
            'mapped_alignment': mapped.get(acc), 'missing_reason': None if acc in mapped else 'no-reported-FHT-alignment'})
    paths = [source, review_path, fht_path, fnsi_path, original_path, coordinates, prior_hits, proteome]
    report = {'schema': 'cannabis-flavone-site-review-v1', 'model_eligible': False,
        'reference_numbering_compatible': True, 'review': review, 'rows': rows,
        'coordinate_command': command, 'coordinate_replay_alignment_count': len(lines),
        'masked_alignments': masked_alignments,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'original_fnsi_leads_retained': len(rows), 'mapped_leads': sum(r['mapped_alignment'] is not None for r in rows),
            'alignments_with_reported_masking': len(masked_alignments), 'new_exact_enzyme_assignments': 0},
        'claim_boundary': 'Local sequence-alignment positional hypotheses only. No vote-based specificity score or transferred enzyme assignment. Reference numbering compatibility does not prove construct identity, structural equivalence, Cannabis product specificity, or in-vivo flux. Original weak/fragment evidence remains in parent annotations.'}
    Path('data/reports/phase1-flavone-site-review.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
