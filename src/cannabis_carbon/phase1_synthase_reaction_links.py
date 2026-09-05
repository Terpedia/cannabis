"""Join revalidated synthase proteins to exact annotated Rhea equations."""
import copy
import hashlib
import json
from pathlib import Path

from .phase1_no_producer_audit import diagnostic_keys
from .phase1_scope import write_rows


def build(network, discovery, search):
    compounds = {c['id']: c for c in network['compounds']}
    reactions = {r['id']: r for r in network['reactions']}
    alignments = {a['id']: a for a in search['passing_alignments']}
    sequences = {p['accession']: p['sequence'] for p in search['cannabis_candidates']}
    references = {p['accession']: p['sequence'] for p in search['reference_sequences']}
    rows, used = [], set()
    for source in discovery['rows']:
        core = reactions[source['reaction_id']]
        for ref in source['reference_matches']:
            accession = ref['accession']
            comments = [c for c in ref['reference_record']['comments'] if c['commentType'] == 'CATALYTIC ACTIVITY']
            for comment in comments:
                for physiology in comment.get('physiologicalReactions', []):
                    crossref = physiology['reactionCrossReference']
                    if crossref['database'] != 'Rhea':
                        continue
                    matches = [(r, s) for r in network['reactions'] for s in r['sources'] if s['source_reaction_id'] == crossref['id']]
                    if len(matches) != 1:
                        raise ValueError('Missing or ambiguous exact Rhea provenance')
                    reaction, provenance = matches[0]
                    if physiology['directionType'] != 'left-to-right':
                        raise ValueError('Unhandled reference physiological orientation')
                    core_source = next(s for s in core['sources'] if s['source_reaction_id'] == source['legacy_source_reaction_id'])
                    comparison = []
                    for role in ('substrates', 'products'):
                        def members(r, p):
                            left = p['source_left_corresponds_to']
                            side = left if role == 'substrates' else ('right' if left == 'left' else 'left')
                            return r[side]
                        a, b = members(core, core_source), members(reaction, provenance)
                        def signature(ms, key=None):
                            return sorted((diagnostic_keys(compounds[m['compound_id']]['smiles'])[key] if key else compounds[m['compound_id']]['smiles'], m['coefficient']) for m in ms)
                        comparison.append({'role': role, 'core_members': copy.deepcopy(a), 'rhea_members': copy.deepcopy(b),
                            'exact_encoded_match': signature(a) == signature(b),
                            'uncharger_only_match': signature(a, 'uncharger') == signature(b, 'uncharger'),
                            'uncharger_and_stereo_removed_match': signature(a, 'uncharger_and_stereo_removed') == signature(b, 'uncharger_and_stereo_removed')})
                    links = []
                    for aid, alignment in sorted(alignments.items()):
                        if alignment['reference_accession'] != accession:
                            continue
                        protein = alignment['cannabis_accession']
                        exact = sequences[protein] == references[accession]
                        links.append({'candidate_accession': protein, 'reference_accession': accession, 'alignment_id': aid,
                            'evidence_class': 'exact-reference-sequence-with-annotated-catalysis' if exact else 'reference-homology-candidate',
                            'sequence_identical': exact, 'direct_candidate_assay_claimed': False})
                    rows.append({'id': accession + ':' + crossref['id'], 'reaction_id': reaction['id'],
                        'core_reaction_id': core['id'], 'annotated_rhea_id': crossref['id'],
                        'canonical_forward_side': provenance['source_left_corresponds_to'],
                        'reference_catalytic_annotation': copy.deepcopy(comment), 'reference_accession': accession,
                        'reference_url': f'https://www.uniprot.org/uniprotkb/{accession}/entry',
                        'source_provenance': copy.deepcopy(provenance), 'core_comparison': comparison,
                        'protein_links': links, 'core_identity_merge_allowed': False})
                    used.update((core['id'], reaction['id']))
    selected = [reactions[rid] for rid in sorted(used)]
    cids = {m['compound_id'] for r in selected for side in ('left', 'right') for m in r[side]}
    return {'schema': 'cannabis-carbon.phase1-synthase-reaction-links.v1', 'rows': rows,
        'reactions': selected, 'compounds': [compounds[c] for c in sorted(cids)],
        'summary': {'exact_rhea_equations': len({r['reaction_id'] for r in rows}),
            'protein_reaction_links': sum(len(r['protein_links']) for r in rows),
            'sequence_identical_links': sum(p['sequence_identical'] for r in rows for p in r['protein_links']),
            'historical_core_equations_merged': 0},
        'claim_boundary': 'Exact directional Rhea provenance joins from reference catalytic annotations. Reference evidence codes and citations are preserved, not independently repeated assays. Sequence identity supports transfer of the reference annotation, not expression, processing, compartment, in-vivo flux or a CO2 route. Homologous proteins remain candidates. Neutralization and stereo-removal comparisons are diagnostics only, not reactions or identity equivalence. No core structures, candidate-model equations or completeness counts changed. Atom tracing remains deferred.'}


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-full-balanced-network', 'phase1-cannabinoid-revalidation-references', 'phase1-cannabinoid-revalidation-search')]
    inputs = [json.loads(p.read_text()) for p in paths]
    if inputs[2]['source_discovery_sha256'] != hashlib.sha256(paths[1].read_bytes()).hexdigest():
        raise ValueError('Discovery lineage changed')
    report = build(*inputs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-synthase-reaction-links.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('link', 'rows'), ('reaction', 'reactions'), ('compound', 'compounds')]
    records = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, collection in groups:
        records.extend((kind, r['id'], r) for r in report[collection])
    count = write_rows(records, sha, Path('data/derived/phase1-synthase-reaction-links.ndjson'))
    print(json.dumps({'summary': report['summary'], 'rows': count, 'sha256': sha}))


if __name__ == '__main__':
    run()
