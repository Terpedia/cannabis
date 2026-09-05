"""Compare FNS-I leads with related parsley enzymes without assigning activity."""
import csv
import hashlib
import json
import subprocess
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from .genome import _fasta
from .phase1_family_search import parse_hits
from .phase1_ureidoglycolate_broad_search import lookup

RAW = Path('data/raw/flavone-fht-comparison')
QUERY = 'ec:1.14.11.9 AND organism_id:4043 AND fragment:false'


def compare(original, hits, reference_ids):
    by_pair = {}
    for group in hits.values():
        for hit in group:
            key = (hit['cannabis_accession'], hit['reference_accession'])
            if key in by_pair:
                raise ValueError('Duplicate comparison alignment')
            by_pair[key] = hit
    rows = []
    for lead in sorted(original['cannabis_candidates'], key=lambda r: r['accession']):
        acc = lead['accession']
        fnsi = [h for h in original['passing_alignments'] if h['cannabis_accession'] == acc]
        comparators = []
        for ref in sorted(reference_ids):
            hit = by_pair.get((acc, ref))
            comparators.append({'reference_accession': ref, 'alignment': hit,
                'status': 'passing' if hit and hit['passes_screen'] else 'weak' if hit else 'no-reported-alignment'})
        rows.append({'accession': acc, 'model_eligible': False,
            'original_fnsi_alignments': fnsi, 'comparators': comparators})
    return rows


def run():
    RAW.mkdir(parents=True, exist_ok=True)
    item = lookup('parsley-fht', QUERY, RAW)
    records = list(csv.DictReader(Path(item['snapshot']).read_text().splitlines(), delimiter='\t'))
    if not records or len({r['Entry'] for r in records}) != len(records):
        raise ValueError('Empty or duplicate comparator inventory')
    annotations = []
    paths = [Path(item['snapshot']), RAW / 'parsley-fht-lookup.json']
    for record in records:
        acc = record['Entry']
        path = RAW / (acc + '.json')
        url = 'https://rest.uniprot.org/uniprotkb/' + acc + '.json'
        if not path.exists():
            with urllib.request.urlopen(url, timeout=45) as response:
                payload = response.read()
            obj = json.loads(payload)
            if obj['primaryAccession'] != acc:
                raise ValueError('Reference identity mismatch')
            path.write_bytes(payload)
        obj = json.loads(path.read_text())
        if obj['primaryAccession'] != acc or obj['organism']['taxonId'] != 4043:
            raise ValueError('Reference identity or organism mismatch')
        annotations.append({'accession': acc, 'query_record': record, 'source_url': url,
            'snapshot': str(path), 'annotation': obj, 'model_eligible': False})
        paths.append(path)
    source = Path('data/reports/phase1-flavone-search.json')
    original = json.loads(source.read_text())
    proteome = Path(original['proteome_path'])
    if hashlib.sha256(proteome.read_bytes()).hexdigest() != original['proteome_sha256']:
        raise ValueError('Changed proteome')
    sequences = _fasta(proteome)
    refs = {a['accession']: a['annotation']['sequence']['value'] for a in annotations}
    fasta = RAW / 'references.fasta'
    fasta.write_text(''.join('>' + acc + '\n' + seq + '\n' for acc, seq in sorted(refs.items())))
    database = RAW / 'references'
    subprocess.run(['diamond', 'makedb', '--in', str(fasta), '--db', str(database)], check=True)
    hits_path = RAW / 'hits.tsv'
    command = list(original['diamond_command'])
    command[command.index('--db') + 1] = str(database)
    command[command.index('--out') + 1] = str(hits_path)
    subprocess.run(command, check=True)
    hits = parse_hits(hits_path.read_text(), sequences.keys(), refs.keys())
    for group in hits.values():
        for hit in group:
            if hit['query_length'] != len(sequences[hit['cannabis_accession']]) or hit['reference_length'] != len(refs[hit['reference_accession']]):
                raise ValueError('Sequence length mismatch')
    rows = compare(original, hits, refs)
    paths += [source, proteome, fasta, hits_path]
    report = {'schema': 'cannabis-flavone-fht-comparison-v1', 'model_eligible': False,
        'generated_at': datetime.now(timezone.utc).isoformat(), 'lookup': item,
        'references': annotations, 'rows': rows, 'all_comparator_alignments': hits,
        'diamond_command': command, 'diamond_version': subprocess.check_output(['diamond', 'version'], text=True).strip(),
        'screen': original['screen'],
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'proteome_sequences': len(sequences), 'comparator_references': len(refs),
            'raw_comparator_alignments': sum(map(len, hits.values())), 'original_fnsi_leads_retained': len(rows),
            'comparator_status_counts': dict(Counter(c['status'] for row in rows for c in row['comparators'])),
            'new_exact_enzyme_assignments': 0},
        'claim_boundary': 'Related-enzyme comparison only. The EC query may recover multifunctional enzymes; original annotations are retained. Neither higher similarity nor absent/weak alignment proves or excludes activity. No comparator is assigned to the missing FNS-I reaction. Original FNS-I evidence and gaps remain unchanged. Study construct numbering must be reconciled before residue transfer.'}
    Path('data/reports/phase1-flavone-fht-comparison.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
