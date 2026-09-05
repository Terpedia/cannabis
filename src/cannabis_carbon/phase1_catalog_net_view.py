"""Publish the canonical catalog-gap report once, as the static graph bundle."""
import hashlib
import json
from pathlib import Path


def run():
    source = Path('data/reports/phase1-catalog-net-gaps.json')
    payload = source.read_bytes(); report = json.loads(payload)
    for path, digest in report['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('Catalog view source checksum mismatch')
    folder = Path('docs/data/catalog-net-view'); folder.mkdir(parents=True, exist_ok=True)
    (folder / 'bundle.json').write_bytes(payload)
    manifest = {'schema': 'cannabis-carbon.phase1-catalog-net-view.v1', 'file': 'bundle.json',
        'bytes': len(payload), 'sha256': hashlib.sha256(payload).hexdigest(),
        'source_sha256': {str(source): hashlib.sha256(payload).hexdigest()}, 'summary': report['summary']}
    (folder / 'index.json').write_text(json.dumps(manifest, separators=(',', ':')) + '\n')
    print(json.dumps(manifest))


if __name__ == '__main__':
    run()
