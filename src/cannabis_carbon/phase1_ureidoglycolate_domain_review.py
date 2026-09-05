"""Audit a weak broad-name lead against explicit reference domain coordinates."""
import hashlib
import json
import subprocess
import urllib.request
from collections import Counter
from pathlib import Path
from .genome import _fasta

RAW = Path('data/raw/ureidoglycolate-domain-review')
QUERY, REFERENCE = 'A0A7J6HP57', 'A0A1I0CMJ9'


def run():
    RAW.mkdir(parents=True, exist_ok=True)
    source = Path('data/reports/phase1-ureidoglycolate-broad-search.json')
    search = json.loads(source.read_text())
    refs = {r['accession']: r['sequence'] for r in search['reference_sequences']}
    sequences = _fasta(Path(search['proteome_path']))
    paths = [source, Path(search['hits_path']), Path(search['proteome_path'])]
    snapshots = {}
    for accession in (QUERY, REFERENCE):
        path = RAW / (accession + '.json')
        if not path.exists():
            with urllib.request.urlopen('https://rest.uniprot.org/uniprotkb/' + accession + '.json', timeout=45) as response:
                payload = response.read()
            data = json.loads(payload)
            if data['primaryAccession'] != accession:
                raise ValueError('Accession mismatch')
            path.write_bytes(payload)
        data = json.loads(path.read_text())
        if data['sequence']['value'] != (sequences[accession] if accession == QUERY else refs[accession]):
            raise ValueError('Annotation sequence differs from searched sequence')
        snapshots[accession] = data
        paths.append(path)
    command = list(search['diamond_command'])
    coords = RAW / 'coordinate-hits.tsv'
    command[command.index('--out') + 1] = str(coords)
    i = command.index('--evalue')
    command[i:i] = ['qstart', 'qend', 'sstart', 'send']
    subprocess.run(command, check=True)
    lines = [line.split('\t') for line in coords.read_text().splitlines()]
    if any(len(f) != 14 for f in lines):
        raise ValueError('Coordinate column mismatch')
    if Counter('\t'.join(f[:10]) for f in lines) != Counter(Path(search['hits_path']).read_text().splitlines()):
        raise ValueError('Coordinate replay differs from original screen')
    selected = [f for f in lines if f[0].split('|')[1] == QUERY and f[1] == REFERENCE]
    if len(selected) != 1:
        raise ValueError('Expected one representative alignment')
    f = selected[0]
    domains = {acc: [x for x in data.get('features', []) if x['type'] == 'Domain'] for acc, data in snapshots.items()}
    overlaps = []
    for acc, begin, end in [(QUERY, int(f[10]), int(f[11])), (REFERENCE, int(f[12]), int(f[13]))]:
        for domain in domains[acc]:
            start, stop = domain['location']['start']['value'], domain['location']['end']['value']
            overlaps.append({'accession': acc, 'domain': domain, 'alignment_start': begin,
                'alignment_end': end, 'overlapping_residues': max(0, min(end, stop) - max(begin, start) + 1)})
    paths.append(coords)
    report = {'schema': 'cannabis-ureidoglycolate-domain-review-v1', 'model_eligible': False,
        'query_accession': QUERY, 'reference_accession': REFERENCE,
        'coordinate_command': command, 'coordinate_replay_alignment_count': len(lines),
        'selected_alignment_columns': f, 'domain_overlaps': overlaps,
        'review_decision': 'Keep the lead ineligible. The selected match overlaps an annotated isopropylmalate-dehydrogenase-like domain in both sequences; the broad reference name is not evidence of ureidoglycolate catalysis by that region.',
        'annotation_snapshots': snapshots,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'claim_boundary': 'One representative weak broad-name hit, not a characterized Cannabis enzyme. Coordinates replay the original full-proteome screen. Domain overlap and automated annotations do not establish ureidoglycolate catalysis or an alternative exact activity. Original thresholds and failed-screen result are unchanged.'}
    Path('data/reports/phase1-ureidoglycolate-domain-review.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps({'alignment': f, 'overlaps': overlaps}), flush=True)


if __name__ == '__main__':
    run()
