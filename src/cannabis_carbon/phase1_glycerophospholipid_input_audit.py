"""Unchanged-model supply audit for all terminal PG/PGP proposal inputs."""
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from .phase1_reaction_completion_net import assemble, validate_certificate
from .phase1_marts_completions import balanced
from .phase1_net_flux import NetModel

BOUNDARY = ('Diagnostic only, not additional metabolite coverage. Individual terminal '
    'precursor supply is tested under the complete unchanged reaction model, inherited '
    'direction exclusions and CO2-only external carbon boundary. Positive results are '
    'exact net-conversion hypotheses with regenerated pre-existing pools allowed. '
    'Individual precursor feasibility is not a joint pathway certificate, and failure '
    'does not prove a unique causal bottleneck or biological absence. Template-branch '
    'input lists are diagnostics, not proof that all alternative routes require them.')


def assemble_current(network, completions, current, layers):
    reactions, compounds, _ = assemble(network, completions)
    for layer in [*layers, current]:
        for c in layer.get('compounds', []):
            if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
                raise ValueError('Compound identity conflict')
            compounds.setdefault(c['id'], c)
        for r in layer.get('added_reactions', []):
            if r['id'] in reactions and any(r[s] != reactions[r['id']][s] for s in ('left','right')):
                raise ValueError('Reaction identity conflict')
            reactions.setdefault(r['id'], r)
    if len(reactions) != current['summary']['balanced_equations']:
        raise ValueError('Incomplete full model reconstruction')
    for r in reactions.values():
        if not balanced([r['left'],r['right']],compounds):
            raise ValueError('Full model balance failure')
    exchange = set(current['external_exchange_compound_ids'])
    if {c for c in exchange if compounds[c]['carbon_count']} != {current['co2_compound_id']}:
        raise ValueError('Changed carbon boundary')
    return reactions, compounds, NetModel(list(reactions.values()),exchange,
        forbidden_step_ids=current['forbidden_step_ids'])


def input_closure(root, candidates):
    producers = defaultdict(list)
    for row in candidates:
        producers[row['compound_id']].append(row)
    pending = [root]; seen = set(); terminal = set(); reaction_ids = set()
    while pending:
        cid = pending.pop()
        if cid in seen:
            continue
        seen.add(cid)
        if not producers[cid]:
            terminal.add(cid)
        for row in producers[cid]:
            reaction_ids.add(row['reaction_id'])
            pending.extend(row['required_precursor_ids'])
    return sorted(terminal), sorted(reaction_ids)


def run():
    root = Path('data/reports')
    current_path = root/'phase1-glycerophospholipid-net.json'
    current = json.loads(current_path.read_bytes())
    paths = [root/'phase1-full-balanced-network.json', root/'phase1-marts-completions.json',
             current_path, root/'phase1-glycerophospholipid-synthesis.json'] + [
                 root/n for n in current['baseline_certificate_reports']]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p,sha in doc.get('source_sha256',{}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed input snapshot')
    network, completions, current, proposals, *layers = docs
    reactions, compounds, model = assemble_current(network,completions,current,layers)
    steps = {s['id']:s for s in model.steps}
    results = []
    for i, row in enumerate(proposals['frontier'],1):
        cid = row['compound_id']
        result = {'compound_id':cid, **model.solve(cid)}
        if result['status'] == 'exact-net-conversion-hypothesis':
            validate_certificate(result,steps,compounds,set(current['external_exchange_compound_ids']),current['co2_compound_id'])
        results.append(result)
        print(f'Input audit {i}/{len(proposals["frontier"])}: {result["status"]}',flush=True)
    result_by_id = {r['compound_id']:r for r in results}
    statuses = {t['cannabisdb_id']:t['net_status'] for t in current['targets']}
    targets = []
    for t in proposals['targets']:
        terminal, rids = input_closure(t['source_form_compound_id'],proposals['precursor_candidates'])
        if not set(terminal) <= result_by_id.keys():
            raise ValueError('Unaudited terminal input')
        targets.append({k:t[k] for k in ('cannabisdb_id','label','compound_id','source_form_compound_id')} | {
            'net_status':statuses[t['cannabisdb_id']], 'terminal_input_ids':terminal,
            'template_reaction_ids':rids, 'unresolved_terminal_input_ids':[c for c in terminal
                if result_by_id[c]['status'] not in ('exact-net-conversion-hypothesis','explicit-exchange-species; not a synthesis target')]})
    used = {s['reaction_id'] for r in results for s in r.get('steps',[])}
    report = {'schema':'cannabis-carbon.phase1-glycerophospholipid-input-audit.v1',
        'claim_boundary':BOUNDARY,'targets':targets,'input_results':results,
        'certificate_reactions':[reactions[r] for r in sorted(used)],
        'compounds':list(compounds.values()),
        'co2_compound_id':current['co2_compound_id'],
        'external_exchange_compound_ids':current['external_exchange_compound_ids'],
        'forbidden_step_ids':current['forbidden_step_ids'],
        'summary':{'full_model_balanced_equations':len(reactions),'target_records':len(targets),
            'terminal_structures':len(results),'terminal_status_counts':dict(Counter(r['status'] for r in results)),
            'infeasible_target_records':sum(t['net_status']=='solver-reported-infeasible' for t in targets),
            'infeasible_targets_with_unresolved_terminal_inputs':sum(t['net_status']=='solver-reported-infeasible' and bool(t['unresolved_terminal_input_ids']) for t in targets),
            'new_CO2_route_claims':0},
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'phase1-glycerophospholipid-input-audit.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps(report['summary']),flush=True)


if __name__ == '__main__':
    run()
