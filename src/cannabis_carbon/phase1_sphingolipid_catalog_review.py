"""Join alternative enzyme hypotheses to exact catalog status, including exclusions."""
import hashlib
import json
from pathlib import Path
from .phase1_reference_discovery import direction_families


def build():
    network_path = Path('data/reports/phase1-full-balanced-network.json')
    search_path = Path('data/reports/phase1-sphingolipid-alternative-search.json')
    discovery_path = Path('data/reports/phase1-sphingolipid-alternative-references.json')
    directions_path = Path('data/raw/phase1-reference-discovery/rhea-directions.tsv')
    network, search, discovery = [json.loads(p.read_text()) for p in (network_path, search_path, discovery_path)]
    if hashlib.sha256(discovery_path.read_bytes()).hexdigest() != search['source_discovery_sha256']:
        raise ValueError('Alternative search lineage mismatch')
    families = direction_families(directions_path.read_text())
    rows = []
    for reference in discovery['rows'][0]['reference_matches']:
        for activity in reference['source_activity']:
            for xref in activity['reaction']['reactionCrossReferences']:
                rid = xref['id']
                if xref['database'] != 'Rhea' or not rid.startswith('RHEA:'):
                    continue
                family = families[rid]
                ids = set(family.values())
                accepted = [r for r in network['reactions'] if any(s['source_reaction_id'].upper() in ids for s in r['sources'])]
                excluded = [r for r in network['excluded_rhea_source_records'] if r['source_reaction_id'].upper() in ids]
                rows.append({'reference_accession': reference['accession'], 'rhea_id': rid,
                    'rhea_family': family, 'source_activity': activity,
                    'passing_alignments': [a for a in search['passing_alignments'] if a['reference_accession'] == reference['accession']],
                    'accepted_balanced_reaction_ids': [r['id'] for r in accepted], 'excluded_catalog_records': excluded,
                    'status': 'catalog-excluded-not-auditable' if excluded and not accepted else 'requires-review',
                    'model_eligible': False})
    return {'schema': 'cannabis-sphingolipid-catalog-review-v1', 'rows': rows,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (network_path, search_path, discovery_path, directions_path)},
        'next_steps': [
            'Resolve specific N-acyl chain and sphingoid-base structures among CannabisDB targets before instantiating any generic ceramide reaction.',
            'Represent every reactant and electron-transfer partner explicitly and audit elements and charge; do not replace cytochrome b5 by an undeclared reducing equivalent.',
            'Assay each Cannabis candidate with defined sphingolipid substrates and distinguish E/Z delta-8 products from delta-6 acyl-CoA products using authentic standards and reference/no-enzyme controls.'
        ],
        'claim_boundary': 'Exact published Rhea-family joins preserve excluded catalog records. Strong reference homology does not turn generic lipid or protein participants into a concrete balanced Cannabis reaction. No target identity or pathway status changed.'}


if __name__ == '__main__':
    Path('data/reports/phase1-sphingolipid-catalog-review.json').write_text(json.dumps(build(), separators=(',', ':')) + '\n')
