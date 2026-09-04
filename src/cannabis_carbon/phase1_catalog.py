"""Deduplicated, full-participant Phase 1 catalog with a complete source ledger."""
import hashlib
import json
from pathlib import Path
from rdkit import RDLogger
from .balance import _reaction_smiles_balance
from .phase1_balance_reference import concrete_participants
from .phase1_overlay import unique_rows


def stable_id(namespace, value):
    return namespace + ':' + hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def assemble(queue, evidence, alternatives, family):
    originals = unique_rows(queue['rows'])
    evidence_rows = unique_rows(evidence['rows'])
    family_rows = unique_rows(family['rows'])
    balanced = {key: row for key, row in originals.items() if row['balance_status'] == 'balanced'}
    if balanced.keys() != evidence_rows.keys() or not family_rows.keys() <= balanced.keys():
        raise ValueError('Evidence scope differs from the balanced reaction set')
    compounds, reactions = {}, {}
    for key, row in balanced.items():
        element, charge = _reaction_smiles_balance(row['reaction_smarts'])
        if not element or not charge or element['status'] != 'balanced' or charge['status'] != 'balanced':
            raise ValueError('A claimed balanced equation fails the current audit')
        participants = concrete_participants(row['reaction_smarts'])
        sides = []
        for side in participants:
            members = []
            for participant in side:
                compound_id = stable_id('structure', participant['smiles'])
                compounds.setdefault(compound_id, {'id': compound_id,
                    **{k: v for k, v in participant.items() if k != 'coefficient'},
                    'identity_basis': 'exact-canonical-isomeric-SMILES; no CannabisDB identity inferred'})
                members.append({'compound_id': compound_id, 'coefficient': participant['coefficient']})
            sides.append(members)
        ev = evidence_rows[key]
        fam = family_rows.get(key, {})
        reactions[key] = {'id': stable_id('reaction-variant', list(key)),
                         'source_reaction_id': key[0], 'reaction_smarts': key[1],
                         'reactants': sides[0], 'products': sides[1],
                         'element_balance': element, 'charge_balance': charge,
                         'equation_orientation': 'as-recorded; not automatically physiological',
                         'required_input_status': 'not-assessed; every reactant is required',
                         'source_urls': row['source_urls'],
                         'source_references': {k: row[k] for k in ('source_uniprot_ids', 'source_genbank_ids', 'source_ec_numbers')},
                         'enzyme_evidence': {
                             'screened_homology_proteins': sorted({h['cannabis_accession'] for h in ev.get('sequence_hits', []) if h['passes_screen']}),
                             'core_reaction_associations': ev.get('core_reaction_evidence', []),
                             'direction_unresolved_family_proteins': fam.get('screened_cannabis_proteins', []),
                             'alignment_evidence_sources': ['phase1-core-enzyme-evidence.json', 'phase1-family-protein-search.json'],
                             'claim_boundary': 'Source references, annotations and homology are separate candidate evidence, not confirmed enzyme activity.'},
                         'alternative_source_links': []}
    alt_by_original = {}
    for row in alternatives['rows']:
        original_key = (row['reaction_id'], row['original_reaction_smarts'])
        if original_key not in originals or original_key in balanced or original_key in alt_by_original:
            raise ValueError('Alternative ledger does not refer uniquely to an excluded source variant')
        links = []
        for candidate in row['balanced_reference_candidates']:
            target = (candidate['rule_id'], candidate['reaction_smarts'])
            if target not in reactions:
                raise ValueError('Alternative is absent from the independently audited balanced catalog')
            link = {'original_reaction_id': original_key[0], 'original_reaction_smarts': original_key[1],
                    'catalog_reaction_id': reactions[target]['id'],
                    'carbon_stoichiometry_matches': candidate['carbon_stoichiometry_matches'],
                    'participant_changes': candidate['participant_changes'],
                    'original_source_urls': originals[original_key]['source_urls'],
                    'original_source_references': {k: originals[original_key][k] for k in ('source_uniprot_ids', 'source_genbank_ids', 'source_ec_numbers')},
                    'association_status': 'candidate-alternative; enzyme evidence not transferred',
                    'required_review': candidate['required_review']}
            reactions[target]['alternative_source_links'].append(link)
            links.append(link)
        alt_by_original[original_key] = links
    if alt_by_original.keys() != originals.keys() - balanced.keys():
        raise ValueError('Alternative ledger must retain every excluded source variant')
    ledger = [{'reaction_id': key[0], 'reaction_smarts': key[1], 'original_balance_status': row['balance_status'],
               'catalog_reaction_id': reactions[key]['id'] if key in reactions else None,
               'balanced_alternative_links': alt_by_original.get(key, []),
               'source_urls': row['source_urls']}
              for key, row in originals.items()]
    return {'schema': 'cannabis-carbon.phase1-reaction-catalog.v1',
            'summary': {'source_variants': len(ledger), 'balanced_reaction_variants': len(reactions),
                        'compound_structures': len(compounds),
                        'source_variants_with_alternative_links': sum(bool(r['balanced_alternative_links']) for r in ledger),
                        'distinct_balanced_reactions_with_alternative_links': sum(bool(r['alternative_source_links']) for r in reactions.values()),
                        'excluded_source_variants': len(originals) - len(balanced),
                        'excluded_source_variants_without_alternatives': sum(not links for links in alt_by_original.values())},
            'compounds': list(compounds.values()), 'reactions': list(reactions.values()), 'source_ledger': ledger,
            'scope': 'Candidate expansion only, not the full Cannabis metabolome or core NetworkDB.',
            'claim_boundary': 'Full equations retain every required reactant and product. Alternative links are not extra reactions or automatic enzyme transfers. Equation orientation and connectivity do not establish physiology or CO2 reachability; no intracellular currencies are seeded.'}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports') / name for name in ('phase1-enzyme-discovery-queue.json',
             'phase1-core-enzyme-evidence.json', 'phase1-balance-reference.json', 'phase1-family-protein-search.json')]
    result = assemble(*(json.loads(path.read_text()) for path in paths))
    result['source_sha256'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    Path('data/reports/phase1-reaction-catalog.json').write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(result['summary'])


def export_table(report_path, output):
    data = report_path.read_bytes()
    report = json.loads(data)
    digest = hashlib.sha256(data).hexdigest()
    records = []
    for kind, collection in [('compound', 'compounds'), ('reaction', 'reactions'), ('source_variant', 'source_ledger')]:
        for row in report[collection]:
            identifier = row.get('id') or stable_id('source-variant', [row['reaction_id'], row['reaction_smarts']])
            records.append({'record_kind': kind, 'record_id': identifier,
                            'source_reaction_id': row.get('source_reaction_id') or row.get('reaction_id'),
                            'record_json': json.dumps(row, separators=(',', ':')), 'report_sha256': digest})
    output.write_text(''.join(json.dumps(row, separators=(',', ':')) + '\n' for row in records))
    return len(records)


if __name__ == '__main__':
    run()
