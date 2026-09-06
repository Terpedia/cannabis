"""Verify public on-demand map assets against one exact Git release manifest."""
import argparse
import hashlib
import json
import subprocess
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath


def safe_path(path):
    p=PurePosixPath(path)
    if not path or p.is_absolute() or '..' in p.parts or str(p)!=path or '?' in path or '#' in path:
        raise ValueError('Unsafe asset path')
    return path


def digest(payload):
    return {'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest()}


def asset_manifest(read_git,folder,extras=()):
    safe_path(folder)
    manifest_path=folder+'/index.json'; manifest_bytes=read_git(manifest_path)
    manifest=json.loads(manifest_bytes)
    if manifest['file']!='bundle.json':
        raise ValueError('Unexpected index file')
    index_path=folder+'/bundle.json'; index_bytes=read_git(index_path)
    if digest(index_bytes)!={k:manifest[k] for k in ('bytes','sha256')}:
        raise ValueError('Committed index manifest mismatch')
    index=json.loads(index_bytes)
    if index['schema']!='cannabis-carbon.sharded-net-view.v1':
        raise ValueError('Wrong view schema')
    result={manifest_path:digest(manifest_bytes),index_path:digest(index_bytes)}
    refs=[index['shared_chemistry'],*index['certificate_files'].values()]
    for ref in refs:
        safe_path(ref['file']); path=folder+'/'+ref['file']
        if path in result:
            raise ValueError('Duplicate asset path')
        if not isinstance(ref['bytes'],int) or ref['bytes']<=0 or len(ref['sha256'])!=64 or any(c not in '0123456789abcdef' for c in ref['sha256']):
            raise ValueError('Invalid asset digest')
        result[path]={k:ref[k] for k in ('bytes','sha256')}
    for path in extras:
        safe_path(path)
        if path in result:
            raise ValueError('Duplicate extra asset')
        result[path]=digest(read_git(path))
    return result


def verify_asset(path,expected,fetch):
    actual=digest(fetch(path,expected['sha256']))
    if actual!=expected:
        raise ValueError('Public asset differs from committed release: '+path)
    return {'path':path,**actual}


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--revision',required=True)
    parser.add_argument('--folder',required=True)
    parser.add_argument('--extra-file',action='append',default=[])
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    revision=subprocess.check_output(['git','rev-parse','--verify',args.revision+'^{commit}'],text=True).strip()
    read_git=lambda p:subprocess.check_output(['git','show',revision+':docs/'+safe_path(p)])
    assets=asset_manifest(read_git,args.folder,args.extra_file)
    origin='https://terpedia.github.io/cannabis/'
    def fetch(path,sha):
        with urllib.request.urlopen(origin+path+'?v='+sha[:16],timeout=55) as response:
            if response.status!=200:
                raise ValueError('Unexpected HTTP status')
            return response.read()
    def check(item):
        return verify_asset(*item,fetch)
    verified=[]
    with ThreadPoolExecutor(max_workers=12) as pool:
        for row in pool.map(check,assets.items()):
            verified.append(row)
            if len(verified)%250==0:
                print(f'Public assets verified: {len(verified)}/{len(assets)}',flush=True)
    report={'schema':'cannabis-carbon.sharded-pages-verification.v1','revision':revision,
        'url':origin,'verified_at':datetime.now(timezone.utc).isoformat(),
        'verified':True,'asset_count':len(verified),'total_bytes':sum(r['bytes'] for r in verified),
        'verification':'Public bytes and SHA256 match the exact committed index manifests for every referenced certificate and shared chemistry asset. Listed extra files and root manifests match direct git-show bytes. This is content verification, not browser visual QA or scientific validation.',
        'assets':verified}
    Path(args.output).write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='assets'}),flush=True)


if __name__=='__main__':
    run()
