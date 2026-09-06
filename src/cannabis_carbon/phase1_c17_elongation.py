"""C15-to-C17 chain-length hypothesis with source-exact carrier/cofactor chemistry."""
import hashlib
import json
from pathlib import Path
from .phase1_odd_chain_elongation import build

BOUNDARY = ('C15-to-C17 substrate-scope hypothesis derived from the Terpedia C18-to-C20 four-step '
    'elongation sequence. Only terminal saturated chain length changes; exact CoA/intermediate '
    'stereochemistry and all cofactors and coefficients are preserved. Cannabis activity and exact '
    'substrate assays are unestablished. Four proposed forward-only steps; no enzyme assignment. '
    'The precursor audit supports a permissive-model net C15 certificate, not physiological supply. '
    'This proposal alone claims no complete pathway or coverage gain.')


def run():
    root = Path('data/reports')
    paths = [root / ('phase1-' + n + '.json') for n in ('odd-chain-net', 'full-balanced-network', 'odd-chain-supply')]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for path, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale source')
    report = build(docs[0], docs[1], shifts=(-3,))
    report['schema'] = 'cannabis-carbon.phase1-c17-elongation.v1'
    report['claim_boundary'] = BOUNDARY
    for r in report['reactions']:
        r['claim_boundary'] = BOUNDARY
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (root / 'phase1-c17-elongation.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
