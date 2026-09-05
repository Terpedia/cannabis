"""Inventory gaps after excluding both unsupported ureide condensations."""
import json
from pathlib import Path
from .phase1_allantoate_alternative_gaps import build


def run():
    prior = json.loads(Path('data/reports/phase1-allantoate-alternative-gaps.json').read_text())
    report = build([p for p in prior['source_sha256'] if p.endswith('search.json')],
                   route_path='data/reports/phase1-ureide-sensitivity.json')
    report['schema'] = 'cannabis-ureide-alternative-gaps-v1'
    Path('data/reports/phase1-ureide-alternative-gaps.json').write_text(
        json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
