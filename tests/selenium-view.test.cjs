const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const {webcrypto,createHash}=require('node:crypto');
const context={crypto:webcrypto,TextDecoder};
vm.runInNewContext(fs.readFileSync('docs/net.js','utf8'),context);
const read=p=>JSON.parse(fs.readFileSync(p));
const current=read('data/reports/phase1-selenium-forward-net.json');
const folder='docs/data/selenium-net-view/';
test('selenium view retains every target and exact historical certificate with validated hashes',()=>{
  const index=read(folder+'bundle.json'), old=read('docs/data/local-speciation-net-view/bundle.json');
  assert.equal(index.targets.length,6220);assert.equal(Object.keys(index.certificate_files).length,2721);
  assert.equal(index.summary.target_status_counts['exact-net-conversion-hypothesis'],2724);
  assert.deepEqual(index.forbidden_step_ids,current.forbidden_step_ids);
  assert.deepEqual(index.targets.map(t=>[t.cannabisdb_id,t.compound_id,t.net_status]),current.targets.map(t=>[t.cannabisdb_id,t.compound_id,t.net_status]));
  for(const [cid,ref] of Object.entries(index.certificate_files)){
    const bytes=fs.readFileSync(folder+ref.file);
    assert.equal(bytes.length,ref.bytes);assert.equal(createHash('sha256').update(bytes).digest('hex'),ref.sha256);
    if(old.certificate_files[cid]) assert.equal(ref.sha256,old.certificate_files[cid].sha256);
    else assert.deepEqual(JSON.parse(bytes),current.new_certificates.find(c=>c.compound_id===cid));
  }
});
test('selenium graph preserves all eleven steps, medium node annotations and forward restriction',async()=>{
  const fetcher=async url=>{const bytes=fs.readFileSync('docs/'+url.split('?')[0]);return {ok:true,json:async()=>JSON.parse(bytes),arrayBuffer:async()=>bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength)};};
  const base=await context.NetView.createLoader(fetcher,'selenium-net-view')();
  const selected=await base.loadTarget('CDB004952');const graph=context.NetView.project(selected,'CDB004952');
  assert.equal(graph.steps.length,11);
  assert.ok(graph.steps.some(s=>s.reaction_id==='selenium-source:KEGG-R03601'&&s.direction_mode==='hypothetical-left-to-right'));
  assert.ok(graph.nodes.some(n=>n.data.compound.smiles==='OO'&&n.data.compound.medium_annotation.review_flags.length));
  assert.ok(graph.edges.every(e=>e.data.required_inputs.length&&e.data.outputs.length));
  const source=graph.steps.find(s=>s.reaction_id==='selenium-source:KEGG-R03601').reaction;
  assert.ok(source.source_name_definition_conflict);assert.equal(source.source_ec,'2.5.1.47');
  const chemistry=read(folder+'chemistry.json');
  const light=chemistry.reactions.filter(r=>r.light_requirement_annotation);
  assert.equal(base.light_reaction_requirements.reactions.length,2);
  assert.deepEqual(read('docs/data/light-reaction-requirements.json'),read('data/curation/light-reaction-requirements.json'));
  const present=base.light_reaction_requirements.reactions.filter(a=>chemistry.reactions.some(r=>r.id===a.model_reaction_id));
  assert.equal(light.length,present.length);
  for(const r of light) assert.equal(r.light_requirement_annotation.energy_inputs[0].entity_type,'photon');
});
