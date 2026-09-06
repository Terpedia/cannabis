import json
import pytest
from cannabis_carbon.verify_sharded_pages import asset_manifest,digest,verify_asset,safe_path


def test_manifest_and_download_validation_fail_closed():
    shared=b'{"reactions":[]}'; cert=b'{"steps":[]}'
    index=json.dumps({'schema':'cannabis-carbon.sharded-net-view.v1',
        'shared_chemistry':{'file':'chemistry.json',**digest(shared)},
        'certificate_files':{'target':{'file':'certificates/a.json',**digest(cert)}}}).encode()
    files={'data/example/bundle.json':index,
        'data/example/index.json':json.dumps({'file':'bundle.json',**digest(index)}).encode(),
        'net.html':b'page'}
    manifest=asset_manifest(files.__getitem__,'data/example',['net.html'])
    assert len(manifest)==5
    row=verify_asset('chemistry.json',digest(shared),lambda p,s:shared)
    assert row['bytes']==len(shared)
    with pytest.raises(ValueError,match='differs'):
        verify_asset('chemistry.json',digest(shared),lambda p,s:cert)
    files['data/example/bundle.json']+=b' '
    with pytest.raises(ValueError,match='manifest mismatch'):
        asset_manifest(files.__getitem__,'data/example')


@pytest.mark.parametrize('path',['../outside','/absolute','a/../b','a?x','a#fragment',''])
def test_unsafe_paths_rejected(path):
    with pytest.raises(ValueError,match='Unsafe'):
        safe_path(path)
