"""Explicit local proton hypotheses, including net-zero zwitterion relocation."""
import hashlib
import json
from collections import Counter
from pathlib import Path

from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced
from .phase1_protonation_audit import fingerprint

BOUNDARY = ('Computed local speciation hypothesis, not a curated proton-transfer reaction '
    'or an exact ChEBI endpoint mapping. Target and reaction-participant identities remain '
    'separate. Independent proton-only normalization preserves heavy atoms, isotopes, '
    'bonds and encoded stereochemistry. Net-zero intramolecular proton relocation is '
    'distinguished from explicit H+ exchange. Reversibility is a sensitivity assumption; '
    'pH, pKa, compartments, catalysis and physiological relevance are unresolved. '
    'A source reaction supports only the partner structure and its equation, not this '
    'speciation proposal. No CO2-route or historical coverage gain is claimed here.')


def proposal(a, b, compounds):
    left, why = fingerprint(compounds[a]['smiles'])
    right, why2 = fingerprint(compounds[b]['smiles'])
    if not left or not right or left['lookup_key'] != right['lookup_key'] or left['canonical_smiles'] == right['canonical_smiles']:
        return None
    delta = right['formal_charge'] - left['formal_charge']
    if right['hydrogen_count'] - left['hydrogen_count'] != delta:
        return None
    sides = [[{'compound_id': a, 'coefficient': 1}], [{'compound_id': b, 'coefficient': 1}]]
    if delta:
        proton = stable_id('structure', '[H+]')
        if proton not in compounds:
            raise ValueError('Missing exact proton')
        sides[0 if delta > 0 else 1].append({'compound_id': proton, 'coefficient': abs(delta)})
    if not balanced(sides, compounds):
        raise ValueError('Speciation equation failed full balance')
    return {'id': stable_id('local-speciation-hypothesis', [a, b]), 'left': sides[0], 'right': sides[1],
        'target_compound_id': a, 'reaction_participant_compound_id': b,
        'hypothesis_type': 'local-speciation',
        'speciation_type': 'explicit-proton-exchange' if delta else 'net-zero-intramolecular-proton-relocation',
        'protons_consumed': delta, 'target_identity_check': left, 'participant_identity_check': right,
        'direction_status': 'explicit-reversible-local-speciation-sensitivity-assumption',
        'balance_status': 'independently-element-isotope-charge-balanced',
        'enzyme_evidence_ids': [], 'claim_boundary': BOUNDARY}


def build(gaps):
    compounds = {c['id']: c for c in gaps['compounds']}
    partners = {r['compound_id']: r for r in gaps['alternative_producers']}
    source_reactions = {r['id']: r for r in gaps['reactions']}
    proposals = {}; excluded = []
    for t in gaps['targets']:
        if t['category'] == 'source-identity-conflict':
            continue
        for cid in t['diagnostic_alternatives']['uncharger']:
            r = proposal(t['compound_id'], cid, compounds)
            if r is None:
                excluded.append({'cannabisdb_id': t['cannabisdb_id'], 'partner_compound_id': cid})
                continue
            rids = sorted({s.rsplit(':', 1)[0] for s in partners[cid]['allowed_producing_step_ids']})
            if not rids or not set(rids) <= source_reactions.keys():
                raise ValueError('Missing partner-producing provenance')
            r.update({'cannabisdb_ids': [], 'source_target_assertions': [],
                'partner_producing_step_ids': partners[cid]['allowed_producing_step_ids'],
                'partner_source_reaction_ids': rids,
                'source_url': 'https://terpedia.github.io/cannabis/data/current-reaction-gaps.json',
                'source_evidence_type': 'computed-local-proton-change-not-curated-reaction'})
            saved = proposals.setdefault(r['id'], r)
            saved['cannabisdb_ids'].append(t['cannabisdb_id'])
            saved['source_target_assertions'].append({k: t[k] for k in ('cannabisdb_id', 'label', 'source_url', 'source_external_ids')})
    used = {p['compound_id'] for r in proposals.values() for side in ('left', 'right') for p in r[side]}
    source_ids = {rid for r in proposals.values() for rid in r['partner_source_reaction_ids']}
    return {'schema': 'cannabis-carbon.phase1-local-speciation-hypotheses.v1',
        'reactions': list(proposals.values()), 'compounds': [compounds[cid] for cid in sorted(used)],
        'partner_source_reactions': [source_reactions[rid] for rid in sorted(source_ids)],
        'excluded_pairs': excluded, 'claim_boundary': BOUNDARY,
        'summary': {'balanced_hypotheses': len(proposals),
            'target_records': len({t for r in proposals.values() for t in r['cannabisdb_ids']}),
            'speciation_type_counts': dict(Counter(r['speciation_type'] for r in proposals.values())),
            'excluded_pairs': len(excluded), 'coverage_gain_claimed': 0}}


def run():
    path = Path('data/reports/phase1-current-reaction-gaps.json')
    gaps = json.loads(path.read_bytes())
    for p, sha in gaps['source_sha256'].items():
        if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
            raise ValueError('Stale source model')
    report = build(gaps)
    report['source_sha256'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()}
    Path('data/reports/phase1-local-speciation-hypotheses.json').write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
