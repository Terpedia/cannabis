"""Build a reaction-level enzyme discovery queue from audited source edges."""
import json
from collections import defaultdict, Counter
from pathlib import Path


def build_queue(expansion_path: Path, balance_path: Path, output: Path):
    expansion = json.loads(expansion_path.read_text())
    audit = json.loads(balance_path.read_text())
    evidence = defaultdict(list)
    for edge in expansion['rows']:
        evidence[(edge.get('reaction_id'), edge.get('reaction_smarts'))].append(edge)
    rows = []
    for reaction in audit['reactions']:
        key = (reaction.get('reaction_id'), reaction.get('reaction_smarts'))
        edges = evidence.get(key, [])
        def values(field):
            return sorted({e[field] for e in edges if e.get(field)})
        uniprot, genbank, ec = values('source_uniprot_id'), values('source_genbank_id'), values('source_ec_number')
        if reaction['status'] != 'balanced':
            action = 'resolve-stoichiometry'
        elif uniprot or genbank:
            action = 'retrieve-reference-sequences-and-search-cannabis-proteome'
        else:
            action = 'find-characterized-reference-enzyme'
        rows.append({
            'reaction_id': key[0], 'reaction_smarts': key[1],
            'balance_status': reaction['status'], 'source_edge_count': len(edges),
            'source_uniprot_ids': uniprot, 'source_genbank_ids': genbank,
            'source_ec_numbers': ec, 'source_urls': values('source_url'),
            'candidate_cannabisdb_ids': sorted({i for e in edges for i in e.get('candidate_cannabisdb_ids', [])}),
            'next_action': action, 'cannabis_enzyme_status': 'unverified',
            'blocker': 'Source enzyme references require species and reaction-specific validation; reference absence means missing evidence.',
        })
    report = {'schema': 'cannabis-carbon.phase1-enzyme-discovery-queue.v1',
              'sources': [str(expansion_path), str(balance_path)],
              'reaction_variant_count': len(rows),
              'action_counts': dict(Counter(r['next_action'] for r in rows)),
              'rows': rows,
              'claim_boundary': 'Each record is a reaction-ID/SMARTS variant. Source protein references are discovery inputs, not established Cannabis enzyme associations. Atom tracing is deferred.'}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    return report['action_counts']


if __name__ == '__main__':
    print(build_queue(Path('data/reports/terpene-identity-set-candidate-expansion.json'),
                      Path('data/reports/terpene-identity-set-candidate-expansion-balance-audit.json'),
                      Path('data/reports/phase1-enzyme-discovery-queue.json')))
