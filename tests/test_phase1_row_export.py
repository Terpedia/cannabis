import unittest
import hashlib
import json
from pathlib import Path
from cannabis_carbon.phase1_row_export import encode, decode


class RowExportTests(unittest.TestCase):
    def test_actual_publication_reconstructs_full_report(self):
        path = Path('data/reports/phase1-triglyceride-net.json')
        rows = [json.loads(line) for line in Path('data/derived/phase1-triglyceride-net.ndjson').read_text().splitlines()]
        self.assertEqual({r['report_sha256'] for r in rows}, {hashlib.sha256(path.read_bytes()).hexdigest()})
        self.assertEqual(decode(list(reversed(rows))), json.loads(path.read_text()))

    def test_lossless_unordered_reassembly(self):
        report = {'schema': 'fixture', 'targets': [{'id': 'same', 'n': 1}, {'id': 'same', 'n': 2}],
                  'empty': [], 'steps': ['forward', 'backward'], 'metadata': {'value': None}}
        self.assertEqual(decode(list(reversed(encode(report, 'a' * 64)))), report)

    def test_missing_duplicate_and_mixed_snapshots_fail(self):
        rows = encode({'targets': [{'id': 'a'}, {'id': 'b'}]}, 'a' * 64)
        with self.assertRaises(ValueError): decode(rows[:-1])
        with self.assertRaises(ValueError): decode(rows + [rows[-1]])
        rows[-1]['report_sha256'] = 'b' * 64
        with self.assertRaises(ValueError): decode(rows)
