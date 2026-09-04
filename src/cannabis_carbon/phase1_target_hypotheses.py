"""Balanced net-production hypotheses; never infer flux from participation."""
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors
from .balance import _reaction_smiles_balance
from .phase1_balance_reference import concrete_participants
from .phase1_catalog import stable_id
from .phase1_core_evidence import rhea_key


def build(coverage, network, expansion):
    core_index = defaultdict(list)
    for reaction in network['reactions']:
        keys = {rhea_key(x) for x in [reaction['id'], *reaction.get('directional_rhea_ids', [])]} - {None}
        for key in keys:
            core_index[key].append(reaction)
    expansion_index = {'expansion:' + r['id']: r for r in expansion['reactions']}
    compounds, reactions, source_index, evidence = {}, {}, {}, {}
    for source in coverage['reaction_ledger']:
        if source['computed_balance_status'] != 'balanced':
            continue
        element, charge = _reaction_smiles_balance(source['reaction_smiles'])
        if not element or not charge or element['status'] != 'balanced' or charge['status'] != 'balanced':
            raise ValueError('Claimed balanced source fails independent audit')
        sides = []
        for participants in concrete_participants(source['reaction_smiles']):
            members = []
            for participant in participants:
                cid = stable_id('structure', participant['smiles'])
                compounds.setdefault(cid, {'id': cid, **{k: v for k, v in participant.items() if k != 'coefficient'}})
                members.append({'compound_id': cid, 'coefficient': participant['coefficient']})
            sides.append(sorted(members, key=lambda m: m['compound_id']))
        flipped = json.dumps(sides[0], sort_keys=True) > json.dumps(sides[1], sort_keys=True)
        if flipped:
            sides.reverse()
        rid = stable_id('balanced-equation', sides)
        reaction = reactions.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
            'balance_status': 'independently-balanced', 'sources': [],
            'direction_status': 'unresolved-in-Cannabis; canonical side ordering is not physiology'})
        support_ids = []
        for core in core_index.get(rhea_key(source['source_reaction_id']), []):
            eid = 'core-evidence:' + core['id']
            evidence.setdefault(eid, {'id': eid, 'core_reaction_id': core['id'],
                'enzyme_ids': core.get('enzyme_ids', []), 'enzyme_associations': core.get('enzyme_associations', []),
                'candidate_proteins': core.get('candidate_proteins', []), 'source_url': core.get('source_url'),
                'source_direction': core.get('direction', {}),
                'join_boundary': 'Explicit source Rhea ID or listed directional-family ID; activity and direction are not automatically transferred.'})
            support_ids.append(eid)
        if source['id'] in expansion_index:
            item = expansion_index[source['id']]
            eid = 'expansion-evidence:' + item['id']
            evidence.setdefault(eid, {'id': eid, 'source_references': item['source_references'],
                'enzyme_evidence': item['enzyme_evidence'], 'source_urls': item['source_urls'],
                'join_boundary': 'Exact expansion source record; homology and family hits remain candidates, not confirmed enzyme activity.'})
            support_ids.append(eid)
        reaction['sources'].append({'coverage_record_id': source['id'],
            'source_reaction_id': source['source_reaction_id'], 'source_layer': source['source_layer'],
            'source_urls': source.get('source_urls') or [source.get('source_url')],
            'source_left_corresponds_to': 'right' if flipped else 'left', 'evidence_ids': support_ids})
        source_index[source['id']] = rid
    # A possible producer needs a positive NET coefficient, not mere right-side presence.
    producers = defaultdict(set)
    net = {}
    for rid, reaction in reactions.items():
        delta = Counter({m['compound_id']: m['coefficient'] for m in reaction['right']})
        delta.subtract({m['compound_id']: m['coefficient'] for m in reaction['left']})
        net[rid] = dict(delta)
        for cid, coefficient in delta.items():
            if coefficient:
                producers[cid].add(rid)  # Explicit all-directions structural scenario only.
    targets, hypotheses = [], []
    for target in coverage['targets']:
        cid = stable_id('structure', target['canonical_isomeric_smiles']) if target['canonical_isomeric_smiles'] else None
        if cid:
            mol = Chem.MolFromSmiles(target['canonical_isomeric_smiles'])
            compounds.setdefault(cid, {'id': cid, 'smiles': target['canonical_isomeric_smiles'],
                'formula': rdMolDescriptors.CalcMolFormula(mol), 'formal_charge': Chem.GetFormalCharge(mol),
                'carbon_count': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms())})
            compounds[cid].setdefault('cannabisdb_ids', []).append(target['cannabisdb_id'])
        matches = {source_index[m['reaction_record_id']] for m in target['reaction_matches']
                   if m['computed_balance_status'] == 'balanced'}
        ids, unchanged = [], []
        for rid in sorted(matches):
            delta = net[rid].get(cid, 0)
            if not delta:
                unchanged.append(rid)
                continue
            reaction = reactions[rid]
            inputs = reaction['left' if delta > 0 else 'right']
            outputs = reaction['right' if delta > 0 else 'left']
            hid = stable_id('net-production-hypothesis', [target['cannabisdb_id'], rid])
            requirements = [{**member, 'availability': 'unestablished',
                'other_structural_producer_count': len(producers[member['compound_id']] - {rid})}
                for member in inputs]
            eid_set = sorted({eid for s in reaction['sources'] for eid in s['evidence_ids']})
            candidate_support = any(has_candidate_support(evidence[eid]) for eid in eid_set)
            blockers = ['physiological-direction-unestablished', 'all-required-inputs-unestablished',
                        'compartment-and-transport-unestablished', 'Cannabis-enzyme-activity-unconfirmed']
            if not candidate_support:
                blockers.append('no-candidate-enzyme-evidence-attached')
            if any(m['compound_id'] == cid for m in inputs):
                blockers.append('requires-target-bootstrap; cannot-establish-de-novo-production')
            if target['structure_status'] == 'stereo-unspecified-or-unknown':
                blockers.append('target-stereochemistry-unresolved')
            hypotheses.append({'id': hid, 'cannabisdb_id': target['cannabisdb_id'], 'compound_id': cid,
                'reaction_id': rid, 'status': 'blocked', 'structural_status': 'balanced-one-step-net-production-hypothesis',
                'direction_mode': 'hypothetical-left-to-right' if delta > 0 else 'hypothetical-right-to-left',
                'net_target_coefficient': abs(delta), 'required_inputs': requirements, 'outputs': outputs,
                'evidence_ids': eid_set, 'has_candidate_enzyme_evidence': candidate_support, 'blockers': blockers,
                'proposed_tests': ['Verify physiological direction and exact substrates/products against primary reaction evidence.',
                    'Assay a Cannabis candidate protein with all listed inputs; confirm target identity against an authentic standard and include no-enzyme controls.' if candidate_support else
                    'Retrieve exact-reaction protein references and screen the Cannabis proteome; then assay candidates with all listed inputs and no-enzyme controls.',
                    'Establish tissue/compartment compatibility and upstream supply for every required input before claiming a pathway.']})
            ids.append(hid)
        status = 'net-production-hypotheses-found' if ids else 'balanced-participation-only-no-net-production' if matches else target['coverage_status']
        targets.append({'cannabisdb_id': target['cannabisdb_id'], 'label': target['label'],
            'compound_id': cid, 'source_smiles': target['source_smiles'], 'structure_status': target['structure_status'],
            'carbon_count': compounds[cid]['carbon_count'] if cid else None,
            'source_url': target['source_url'], 'status': status, 'hypothesis_ids': ids,
            'unchanged_participant_reaction_ids': unchanged, 'coverage_status': target['coverage_status'],
            'next_step': 'Review inorganic nutrient uptake, transport and required-input supply; this record contains no carbon.' if cid and compounds[cid]['carbon_count'] == 0 else
            'Resolve direction, enzyme and input-supply blockers.' if ids else
            'Curate an exact-identity balanced producing reaction; absence from this catalog is not biological absence.'})
    return {'schema': 'cannabis-carbon.phase1-target-hypotheses.v1',
        'summary': {'cannabisdb_records': len(targets), 'target_status_counts': dict(Counter(t['status'] for t in targets)),
            'exact_compound_structures': len(compounds), 'deduplicated_balanced_equations': len(reactions),
            'balanced_source_records': len(source_index), 'one_step_hypotheses': len(hypotheses),
            'carbon_bearing_target_records': sum(bool(t['carbon_count']) for t in targets),
            'carbon_bearing_target_status_counts': dict(Counter(t['status'] for t in targets if t['carbon_count'])),
            'carbon_free_target_records': sum(t['carbon_count'] == 0 for t in targets),
            'carbon_bearing_target_hypotheses': sum(compounds[h['compound_id']]['carbon_count'] > 0 for h in hypotheses),
            'hypotheses_with_candidate_enzyme_evidence': sum(h['has_candidate_enzyme_evidence'] for h in hypotheses),
            'hypotheses_without_candidate_enzyme_evidence': sum(not h['has_candidate_enzyme_evidence'] for h in hypotheses),
            'hypotheses_requiring_target_bootstrap': sum(any('requires-target-bootstrap' in b for b in h['blockers']) for h in hypotheses)},
        'targets': targets, 'compounds': list(compounds.values()), 'reactions': list(reactions.values()),
        'hypotheses': hypotheses, 'enzyme_evidence': list(evidence.values()),
        'scope': 'Every CannabisDB record retained. Balanced equations from the full-Rhea target-participation audit only; exact full-equation duplicates including reverse encodings are merged, with every source orientation retained.',
        'claim_boundary': 'One-step hypotheses are not complete pathways. Both directions are examined only as a structural scenario, not physiological reversibility. Every input is required; no intracellular currencies or other seeds are supplied. Presence on both sides with zero net production gives no producing hypothesis. Annotation, source reference and homology are not confirmed Cannabis activity. Atom tracing is deferred.'}


def has_candidate_support(record):
    ev = record.get('enzyme_evidence', {})
    return bool(record.get('enzyme_ids') or record.get('enzyme_associations') or record.get('candidate_proteins')
                or ev.get('screened_homology_proteins') or ev.get('direction_unresolved_family_proteins')
                or any(r.get('enzyme_association_ids') or r.get('candidate_proteins')
                       for r in ev.get('core_reaction_associations', [])))


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-target-rhea-coverage.json'), Path('docs/data/networkdb.json'),
             Path('data/reports/phase1-reaction-catalog.json')]
    sources = [json.loads(p.read_text()) for p in paths]
    for path in paths[1:]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != sources[0]['source_sha256'][str(path)]:
            raise ValueError('Coverage source checksum mismatch')
    result = build(*sources)
    result['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result['rdkit_version'] = rdBase.rdkitVersion
    output = Path('data/reports/phase1-target-hypotheses.json')
    output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(json.dumps(result['summary']))


def export_table(report_path, output):
    raw = report_path.read_bytes()
    report = json.loads(raw)
    digest = hashlib.sha256(raw).hexdigest()
    count = 0
    with output.open('w') as stream:
        for kind, collection in [('target', 'targets'), ('compound', 'compounds'), ('reaction', 'reactions'),
                                 ('hypothesis', 'hypotheses'), ('enzyme_evidence', 'enzyme_evidence')]:
            for row in report[collection]:
                stream.write(json.dumps({'record_kind': kind, 'record_id': row.get('id') or row['cannabisdb_id'],
                    'status': row.get('status'), 'record_json': json.dumps(row, separators=(',', ':')),
                    'report_sha256': digest}, separators=(',', ':')) + '\n')
                count += 1
    return count


if __name__ == '__main__':
    run()
