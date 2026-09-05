"""Separate exact reaction gaps from diagnostic identity leads, without merging."""
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

from .phase1_scope import write_rows
from .phase1_target_coverage import encoded_structure


def diagnostic_keys(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or any(a.GetAtomicNum() == 0 for a in molecule.GetAtoms()):
        raise ValueError('Concrete target structure required')
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    def key(mol):
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    no_stereo = Chem.Mol(molecule)
    Chem.RemoveStereochemistry(no_stereo)
    uncharged = rdMolStandardize.Uncharger().uncharge(molecule)
    both = Chem.Mol(uncharged)
    Chem.RemoveStereochemistry(both)
    return {'stereo_removed': key(no_stereo), 'uncharger': key(uncharged),
            'uncharger_and_stereo_removed': key(both)}


def build(network, candidates, source_compounds, xml_table):
    compounds = {c['id']: c for c in network['compounds']}
    source = {c['id']: c for c in source_compounds['compounds']}
    xml = {r['accession']: r for r in xml_table['rows']}
    if len(xml) != len(xml_table['rows']) or set(xml) != set(source):
        raise ValueError('XML accession inventory mismatch')
    conflicts = []
    for accession, sdf in source.items():
        counterpart = xml[accession]
        if encoded_structure(sdf['smiles'])[0] != encoded_structure(counterpart['smiles'])[0]:
            conflicts.append({'cannabisdb_id': accession,
                'sdf_derived_assertion': {k: sdf.get(k) for k in ('label', 'smiles', 'formula', 'inchikey', 'external_ids')},
                'xml_assertion': {k: counterpart.get(k) for k in ('name', 'smiles', 'formula', 'inchi', 'inchikey', 'iupac_name', 'external_ids')},
                'status': 'source-structure-disagreement; identity-resolution-required',
                'claim_boundary': 'Neither source is automatically selected as correct. XML nomenclature and external IDs may themselves disagree with its structure fields. Preserve assertions and resolve against primary structure evidence before changing carbon accounting.'})
    conflict_ids = {r['cannabisdb_id'] for r in conflicts}
    scenarios = candidates['scenarios']
    gaps = [t for t in scenarios[0]['targets'] if t['net_status'] == 'no-net-producing-candidate-equation']
    gap_ids = {t['cannabisdb_id'] for t in gaps}
    if any({t['cannabisdb_id'] for t in s['targets'] if t['net_status'] == 'no-net-producing-candidate-equation'} != gap_ids for s in scenarios):
        raise ValueError('No-producer inventory differs between scenarios')
    if len(gaps) != len(gap_ids) or set(source) != {t['cannabisdb_id'] for t in network['targets']}:
        raise ValueError('Duplicate or incomplete target inventory')
    participation, producers = defaultdict(set), defaultdict(set)
    reaction_by_id = {r['id']: r for r in network['reactions']}
    for reaction in network['reactions']:
        net = defaultdict(Fraction)
        for side, sign in [('left', -1), ('right', 1)]:
            for member in reaction[side]:
                cid = member['compound_id']
                participation[cid].add(reaction['id'])
                net[cid] += sign * Fraction(member['coefficient'])
        # Both orientations are considered only for this catalog-level upper bound.
        for cid, amount in net.items():
            if amount:
                producers[cid].add(reaction['id'])
    indexes = {kind: defaultdict(set) for kind in ('stereo_removed', 'uncharger', 'uncharger_and_stereo_removed')}
    for cid in sorted(participation):
        for kind, key in diagnostic_keys(compounds[cid]['smiles']).items():
            indexes[kind][key].add(cid)
    rows, used_alternatives, used_reactions = [], set(), set()
    for target in gaps:
        cid = target['compound_id']
        exact_producers = sorted(producers[cid])
        if set(exact_producers) & candidates['candidate_reaction_evidence_ids'].keys():
            raise ValueError('No-producer target has an admitted producer')
        # Diagnostics never overwrite the source identity or create graph edges.
        alternatives = {kind: sorted(indexes[kind].get(key, set()) - {cid})
                        for kind, key in diagnostic_keys(compounds[cid]['smiles']).items()}
        status = ('exact-catalog-producer-without-candidate' if exact_producers else
                  'exact-catalog-participant-with-zero-net-only' if participation[cid] else
                  'diagnostic-identity-leads-only' if any(alternatives.values()) else
                  'no-exact-or-diagnostic-catalog-match')
        original = source[target['cannabisdb_id']]
        row = {k: target[k] for k in ('cannabisdb_id', 'compound_id', 'label', 'carbon_count')}
        row.update(status=status, source_smiles=original['smiles'], canonical_smiles=compounds[cid]['smiles'],
                   source_url=original.get('source_url'), source_external_ids=original.get('external_ids', {}),
                   exact_catalog_participation_reaction_ids=sorted(participation[cid]),
                   exact_catalog_net_producer_reaction_ids=exact_producers, diagnostic_alternatives=alternatives,
                   source_identity_status='source-structure-disagreement' if target['cannabisdb_id'] in conflict_ids else 'SDF-derived-and-XML-encoded-structures-agree',
                   next_action=('Resolve SDF/XML structure disagreement before pathway inference or carbon-accounting changes.' if target['cannabisdb_id'] in conflict_ids else
                                'Review exact reaction enzyme evidence and all-input supply.' if exact_producers else
                                'Review identity alternatives against source records; do not merge.' if any(alternatives.values()) else
                                'Find or infer a complete balanced producing reaction with explicit provenance.'))
        rows.append(row)
        used_reactions.update(participation[cid])
        used_alternatives.update(a for matches in alternatives.values() for a in matches)
    for cid in used_alternatives:
        used_reactions.update(participation[cid])
    return {'schema': 'cannabis-carbon.phase1-no-producer-audit.v1', 'targets': rows,
            'source_identity_conflicts': conflicts,
            'alternative_compounds': [{**compounds[cid], 'catalog_participation_reaction_ids': sorted(participation[cid]),
                                      'catalog_net_producer_reaction_ids': sorted(producers[cid])} for cid in sorted(used_alternatives)],
            'reactions': [reaction_by_id[rid] for rid in sorted(used_reactions)],
            'summary': {'all_target_records': len(network['targets']), 'no_candidate_producer_records': len(rows),
                        'exact_target_structures': len({r['compound_id'] for r in rows}),
                        'status_counts': dict(Counter(r['status'] for r in rows)),
                        'diagnostic_target_counts': {kind: sum(bool(r['diagnostic_alternatives'][kind]) for r in rows) for kind in indexes},
                        'alternative_compound_structures': len(used_alternatives)},
            'rdkit_version': rdBase.rdkitVersion,
            'source_identity_summary': {'compared_accessions': len(source), 'structure_disagreement_records': len(conflicts),
                                        'no_producer_records_with_structure_disagreement': len(gap_ids & conflict_ids)},
            'claim_boundary': 'All no-producer target records retained, based on the current SDF-derived structures, including explicitly flagged XML disagreements. Exact catalog producers are chemistry upper bounds with unresolved Cannabis evidence and direction, not pathways. Stereo-removal and Uncharger fingerprints supply diagnostic alternatives only: these transformations are not identity equivalence, proton-transfer reactions, enzyme assignments or permission to merge stereoisomers, charges, isotopes, salts or tautomers. Isotope labels are retained in diagnostic keys. No target, reaction, pathway or completeness count is promoted. Source xrefs are retained as source claims, not validated identity joins. Atom tracing remains deferred.'}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path(p) for p in ('data/reports/phase1-full-balanced-network.json',
             'data/reports/phase1-replacement-candidate-net.json', 'docs/data/compounds.json',
             'data/terpedia/cannabisdb-compounds.json')]
    inputs = [json.loads(p.read_text()) for p in paths]
    for report in inputs[:2]:
        for path, sha in report.get('source_sha256', {}).items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                raise ValueError('Source lineage changed')
    report = build(*inputs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    xml_snapshot = Path('data/terpedia/cannabisdb-compounds.xml.gz')
    report['source_sha256'][str(xml_snapshot)] = hashlib.sha256(xml_snapshot.read_bytes()).hexdigest()
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-no-producer-audit.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('target', 'targets', 'cannabisdb_id'), ('alternative_compound', 'alternative_compounds', 'id'), ('reaction', 'reactions', 'id'), ('source_identity_conflict', 'source_identity_conflicts', 'cannabisdb_id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, key, identifier in groups:
        rows.extend((kind, r[identifier], r) for r in report[key])
    count = write_rows(rows, sha, Path('data/derived/phase1-no-producer-audit.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
