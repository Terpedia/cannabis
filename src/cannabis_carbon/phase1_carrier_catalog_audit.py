"""Catalog-wide carrier-source joins; no repaired pathway or balance claim."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from .phase1_glycerophospholipid_input_audit import assemble_current

RH = 'http://rdf.rhea-db.org/'


def index_sources(properties, occurrences):
    triples = defaultdict(list)
    for raw in properties:
        triples[raw['subject']['value']].append({'predicate': raw['predicate']['value'], 'object': raw['object']})
    carrier_ids = {s for s, rows in triples.items() if any(r['predicate'] == RH + 'reactivePart' for r in rows)}
    carriers = []
    for cid in sorted(carrier_ids):
        parts = [r['object']['value'] for r in triples[cid] if r['predicate'] == RH + 'reactivePart']
        if any(p not in triples for p in parts):
            raise ValueError('Missing reactive-part definition')
        carriers.append({'id': cid, 'properties': triples[cid],
            'reactive_parts': [{'id': p, 'properties': triples[p]} for p in parts]})
    by_reaction = defaultdict(list)
    for raw in occurrences:
        row = {k: v['value'] for k, v in raw.items()}
        if row['carrier'] not in carrier_ids:
            raise ValueError('Unknown source carrier')
        rid = 'RHEA:' + row['reaction'].removeprefix(RH)
        by_reaction[rid].append(row)
    return carriers, by_reaction


def context_key(rows):
    # Ignore orientation only for detecting context collapse, not reaction execution.
    return tuple(sorted(tuple(sorted((r['carrier'], r['coefficientPredicate']) for r in rows
        if r['sideRole'] == RH + role)) for role in ('substrates', 'products')))


def run():
    root = Path('data/reports'); raw = Path('data/raw/rhea-carrier-snapshot-20260906')
    current_path = root / 'phase1-geranial-reduction-net.json'
    current = json.loads(current_path.read_bytes())
    paths = [current_path, root / 'phase1-full-balanced-network.json', root / 'phase1-marts-completions.json',
        raw / 'manifest.json', raw / 'carrier-properties.json', raw / 'carrier-occurrences.json',
        raw / 'independent-counts.json'] + [root / n for n in current['baseline_certificate_reports']]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale model source')
    for receipt in docs[3]['retrievals']:
        if hashlib.sha256(Path(receipt['path']).read_bytes()).hexdigest() != receipt['sha256']:
            raise ValueError('Stale carrier snapshot')
    carriers, source_index = index_sources(docs[4]['results']['bindings'], docs[5]['results']['bindings'])
    expected = {r['name']: r['result'] for r in docs[6]['counts']}
    if len(carriers) != expected['distinct-carriers-with-reactive-parts'] or sum(map(len, source_index.values())) != expected['directed-carrier-occurrence-rows']:
        raise ValueError('Incomplete source snapshot')
    reactions, _, _ = assemble_current(docs[1], docs[2], current, docs[7:])
    equations = []
    for rid, reaction in reactions.items():
        joins = []
        for source in reaction.get('sources', []):
            sid = source.get('source_reaction_id')
            if sid in source_index:
                joins.append({'source_record': source, 'carrier_occurrences': source_index[sid]})
        if joins:
            signatures = {context_key(j['carrier_occurrences']) for j in joins}
            equations.append({'id': rid, 'model_reaction': reaction, 'carrier_source_joins': joins,
                'distinct_unoriented_carrier_contexts': len(signatures),
                'requires_carrier_specific_reconstruction': True})
    flagged = {r['id'] for r in equations}
    certs = {c['compound_id']: c for doc in [*docs[7:], current]
        for c in doc.get('certificates', []) + doc.get('new_certificates', [])}
    if set(certs) != {t['compound_id'] for t in current['targets'] if t['net_status'] == 'exact-net-conversion-hypothesis'}:
        raise ValueError('Incomplete certificate set')
    hits = {cid: sorted({s['reaction_id'] for s in c['steps']} & flagged) for cid, c in certs.items()}
    targets = [{k: t[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'net_status')} | {
        'carrier_source_equation_ids': hits.get(t['compound_id'], []),
        'carrier_audit_status': 'carrier-source-reconstruction-required' if hits.get(t['compound_id']) else
            'no-direct-carrier-source-hit-not-full-identity-validation' if t['compound_id'] in certs else 'no-saved-certificate-to-screen'}
        for t in current['targets']]
    report = {'schema': 'cannabis-carbon.phase1-carrier-catalog-audit.v1', 'carriers': carriers,
        'equations': equations, 'targets': targets,
        'summary': {'target_records': len(targets), 'model_equations': len(reactions),
            'source_carriers': len(carriers), 'source_directed_reactions_with_carriers': len(source_index),
            'source_occurrences': sum(map(len, source_index.values())),
            'model_equations_with_direct_carrier_source': len(equations),
            'model_equations_merging_multiple_carrier_contexts': sum(r['distinct_unoriented_carrier_contexts'] > 1 for r in equations),
            'screened_certificates': len(certs), 'affected_certificates': sum(bool(v) for v in hits.values()),
            'affected_inventory_records': sum(bool(t['carrier_source_equation_ids']) for t in targets),
            'corrected_pathway_coverage_established': False},
        'claim_boundary': 'Direct-source join across the complete current model and all inventory records. '
            'Current Rhea metadata and historical SMILES snapshots may differ in release. Source carrier presence requires '
            'identity-aware reconstruction; balance of reactive-part SMILES does not validate the full participants. '
            'No-hit results are not clearance: polymer-only contexts, inferred completions, template-derived hypotheses and '
            'sources without retained Rhea IDs need separate provenance audits. No reaction is repaired or removed here; '
            'no new full-pathway, physiological-direction, enzyme, nutrient or growth claim is made.',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root / 'phase1-carrier-catalog-audit.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
