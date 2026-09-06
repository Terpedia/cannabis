const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const {webcrypto,createHash}=require('node:crypto');
const context={crypto:webcrypto,TextDecoder};
vm.runInNewContext(fs.readFileSync('docs/net.js','utf8'),context);
const read=p=>JSON.parse(fs.readFileSync(p));
const report=read('data/reports/phase1-ketone-stereo-net.json');
const folder='docs/data/ketone-stereo-net-view/';
test('ketone map preserves all historical certificates and full target inventory',()=>{
  const index=read(folder+'bundle.json'), old=read('docs/data/selenium-net-view/bundle.json');
  assert.equal(index.targets.length,6220);
  assert.equal(Object.keys(index.certificate_files).length,2727);
  assert.equal(index.summary.target_status_counts['exact-net-conversion-hypothesis'],2730);
  assert.deepEqual(index.forbidden_step_ids,report.forbidden_step_ids);
  assert.deepEqual(index.targets.map(t=>[t.cannabisdb_id,t.compound_id,t.net_status]),report.targets.map(t=>[t.cannabisdb_id,t.compound_id,t.net_status]));
  assert.match(index.view_boundary,/reaction-class analogies/);
  assert.match(index.medium_annotation_scope,/not recomputed/);
  for(const [cid,ref] of Object.entries(index.certificate_files)){
    const bytes=fs.readFileSync(folder+ref.file);
    assert.equal(bytes.length,ref.bytes);
    assert.equal(createHash('sha256').update(bytes).digest('hex'),ref.sha256);
    if(old.certificate_files[cid]) assert.equal(ref.sha256,old.certificate_files[cid].sha256);
    else assert.deepEqual(JSON.parse(bytes),report.new_certificates.find(c=>c.compound_id===cid));
  }
});
test('all six new graphs retain explicit redox hypotheses, inputs and direction constraints',async()=>{
  const fetcher=async url=>{const bytes=fs.readFileSync('docs/'+url.split('?')[0]);return {ok:true,json:async()=>JSON.parse(bytes),arrayBuffer:async()=>bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength)};};
  const base=await context.NetView.createLoader(fetcher,'ketone-stereo-net-view')();
  for(const t of report.targets.filter(t=>t.new_certificate)){
    const selected=await base.loadTarget(t.cannabisdb_id);
    const graph=context.NetView.project(selected,t.cannabisdb_id);
    const cert=report.new_certificates.find(c=>c.compound_id===t.compound_id);
    assert.equal(graph.steps.length,cert.steps.length);
    const added=graph.steps.filter(s=>s.reaction_id.startsWith('ketone-stereo-hypothesis:'));
    assert.ok(added.length);
    for(const s of added){
      assert.equal(s.direction_mode,'hypothetical-left-to-right');
      assert.equal(s.reaction.enzyme_evidence_ids.length,0);
      assert.equal(s.reaction.is_route_sensitivity,true);
      assert.match(s.reaction.source_evidence_type,/reaction-class-analogy/);
      assert.match(s.reaction.claim_boundary,/not substrate scope/);
    }
    assert.ok(graph.edges.every(e=>e.data.required_inputs.length&&e.data.outputs.length));
  }
});
