"""Review-only exact PG/PGP/CDP-DAG charge-state gaps, without stereo filling."""
import hashlib
import json
from pathlib import Path
from rdkit import RDLogger
from .phase1_amino_phospholipid_speciation import build as build_speciation

RULES = {'PG': 'RHEA:33752', 'PGP': 'RHEA:12594', 'CDP-DAG': 'RHEA:16230'}
BOUNDARY = ('Review-only generic-source-scaffold proton-exchange hypothesis, not an '
    'exact ChEBI mapping or curated protonation reaction. Neutral and source-charge '
    'structures remain separate identities. Heavy atoms, bonds, isotopes and encoded '
    'stereochemistry are unchanged; unspecified stereochemistry is never filled. '
    'Source reaction provenance establishes the generic product scaffold only, not '
    'exact lipid specificity, pKa, tissue pH, Cannabis enzyme activity, precursor '
    'supply or CO2 conversion. This audit adds no pathways to the model.')


def build(parent, catalog):
    return build_speciation(parent, catalog, rules=RULES, boundary=BOUNDARY,
        schema='cannabis-carbon.phase1-glycerophospholipid-speciation.v1')


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-amino-phospholipid-net.json'),
             Path('data/raw/phase1-balance-reference-catalog.json')]
    report = build(*(json.loads(p.read_bytes()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-glycerophospholipid-speciation.json').write_text(
        json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
