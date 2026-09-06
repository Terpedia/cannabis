"""Static certificate-on-demand packaging without omitting equation participants."""
import hashlib
import json
from pathlib import Path


def write_json(path, value):
    payload = (json.dumps(value, separators=(',', ':'))+'\n').encode()
    path.write_bytes(payload)
    return {'file':path.name, 'sha256':hashlib.sha256(payload).hexdigest(), 'bytes':len(payload)}


def write_view(bundle, folder):
    """Shared chemistry stays intact; only the selected certificate is downloaded."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    certificates = {c['compound_id']:c for c in bundle['certificates']}
    reactions = {r['id']:r for r in bundle['reactions']}
    compounds = {c['id']:c for c in bundle['compounds']}
    evidence = {e['id']:e for e in bundle['enzyme_evidence']}
    if len(certificates)!=len(bundle['certificates']) or len(reactions)!=len(bundle['reactions']):
        raise ValueError('Duplicate certificate or reaction identity')
    if len(compounds)!=len(bundle['compounds']) or len(evidence)!=len(bundle['enzyme_evidence']):
        raise ValueError('Duplicate compound or evidence identity')
    if len({t['cannabisdb_id'] for t in bundle['targets']})!=len(bundle['targets']):
        raise ValueError('Duplicate target accession')
    for r in reactions.values():
        if not set(r['enzyme_evidence_ids']) <= evidence.keys():
            raise ValueError('Missing reaction evidence')
        if not {p['compound_id'] for side in ('left','right') for p in r[side]} <= compounds.keys():
            raise ValueError('Missing full-equation participant')
    gaps = {}
    for cid, cert in certificates.items():
        rids = {s['reaction_id'] for s in cert['steps']}
        if not rids <= reactions.keys():
            raise ValueError('Missing certificate reaction')
        gaps[cid] = sorted(rid for rid in rids if not reactions[rid]['enzyme_evidence_ids'])
    targets = []
    for t in bundle['targets']:
        cid = t.get('certificate_compound_id')
        if cid is not None and (cid not in certificates or cid != t['compound_id']):
            raise ValueError('Target certificate identity mismatch')
        actual = gaps.get(cid, [])
        if t.get('missing_candidate_reaction_ids', []) != actual:
            raise ValueError('Target enzyme-gap list mismatch')
        targets.append({k:v for k,v in t.items() if k!='missing_candidate_reaction_ids'} |
                       {'missing_candidate_reaction_count':len(actual)})
    if {t['certificate_compound_id'] for t in targets if t.get('certificate_compound_id')} != certificates.keys():
        raise ValueError('Orphan certificate')
    shared = {k:bundle[k] for k in ('reactions','compounds','enzyme_evidence')}
    shared_ref = write_json(folder/'chemistry.json', shared)
    routes = {}
    route_folder = folder/'certificates'
    route_folder.mkdir(exist_ok=True)
    for cid, cert in certificates.items():
        # Content-independent identity filenames are versioned by full content hash.
        name = hashlib.sha256(cid.encode()).hexdigest()+'.json'
        ref = write_json(route_folder/name, cert)
        routes[cid] = {**ref, 'file':'certificates/'+name}
    index = {k:v for k,v in bundle.items() if k not in ('targets','certificates','reactions','compounds','enzyme_evidence')}
    index.update({'schema':'cannabis-carbon.sharded-net-view.v1', 'source_bundle_schema':bundle['schema'],
                  'targets':targets, 'certificate_files':routes, 'shared_chemistry':shared_ref,
                  'loading_boundary':'All target records are indexed. Full chemistry is shared; exact certificates load on demand. No equation participants or evidence are discarded.'})
    manifest = write_json(folder/'bundle.json', index)
    manifest['schema'] = index['schema']
    write_json(folder/'index.json', manifest)
    return {'target_records':len(targets), 'certificates':len(routes),
            'index_bytes':manifest['bytes'], 'shared_chemistry_bytes':shared_ref['bytes'],
            'largest_certificate_bytes':max((r['bytes'] for r in routes.values()),default=0),
            'total_certificate_bytes':sum(r['bytes'] for r in routes.values())}
