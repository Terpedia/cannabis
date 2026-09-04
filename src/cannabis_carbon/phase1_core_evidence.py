"""Attach explicit Rhea-family core associations to Phase 1 search gaps."""
import json
import hashlib
from collections import defaultdict
from pathlib import Path


def rhea_key(value):
    value = str(value)
    if value.lower().startswith('rhea:'):
        value = value[5:]
    return value if value.isdigit() else None


def attach(search_path, network_path, output):
    search = json.loads(search_path.read_text())
    network = json.loads(network_path.read_text())
    index = defaultdict(list)
    for reaction in network['reactions']:
        keys = {rhea_key(x) for x in [reaction['id'], *(reaction.get('directional_rhea_ids') or [])]}
        for key in keys - {None}:
            index[key].append(reaction)
    rows = []
    for row in search['rows']:
        matches = index.get(rhea_key(row['reaction_id']), [])
        associations = []
        for reaction in matches:
            associations.append({
                'core_reaction_id': reaction['id'],
                'matched_source_reaction_id': row['reaction_id'],
                'join_method': 'explicit-core-reaction-or-listed-directional-Rhea-ID',
                'source_url': reaction.get('source_url'),
                'equation': reaction.get('equation'), 'direction': reaction.get('direction'),
                'enzyme_association_ids': reaction.get('enzyme_ids', []),
                'enzyme_associations': reaction.get('enzyme_associations', []),
                'reactants': reaction.get('reactants', []),
                'products': reaction.get('products', []),
                'participant_orientation': 'core-equation-orientation-not-expansion-direction',
                'candidate_proteins': reaction.get('candidate_proteins', []),
                'ec_numbers': reaction.get('ec_numbers', []),
                'claim_boundary': 'Core enzyme associations may be annotation-derived. Rhea-family membership does not establish physiological direction or activity for a stereo-insensitive expansion endpoint.',
            })
        rows.append({**row, 'core_reaction_evidence': associations})
    summary = {
        'balanced_variants': len(rows),
        'variants_with_core_reaction_evidence': sum(bool(r['core_reaction_evidence']) for r in rows),
        'missing_reference_variants_with_core_evidence': sum(r['search_status'] == 'no-reference-sequence' and bool(r['core_reaction_evidence']) for r in rows),
        'variants_with_core_enzyme_associations': sum(any(a['enzyme_association_ids'] for a in r['core_reaction_evidence']) for r in rows),
    }
    output.write_text(json.dumps({'schema': 'cannabis-carbon.phase1-core-enzyme-evidence.v1',
        'source_search': str(search_path), 'source_network': str(network_path),
        'source_sha256': {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (search_path, network_path)},
        'summary': summary, 'rows': rows,
        'claim_boundary': 'Core association evidence supplements the sequence search. Original search statuses remain unchanged; associations are not biochemical confirmation.'}, separators=(',', ':')) + '\n')
    return summary


if __name__ == '__main__':
    print(attach(Path('data/reports/phase1-targeted-protein-search.json'),
                 Path('docs/data/networkdb.json'), Path('data/reports/phase1-core-enzyme-evidence.json')))
