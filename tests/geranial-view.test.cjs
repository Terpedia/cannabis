const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const {webcrypto,createHash}=require('node:crypto');
const read=p=>JSON.parse(fs.readFileSync(p));
const report=read('data/reports/phase1-geranial-reduction-net.json');
const audit=read('data/reports/phase1-geranial-evidence-audit.json');
const folder='docs/data/geranial-net-view/';
const fetcher=async url=>{const b=fs.readFileSync('docs/'+url.split('?')[0]);return {ok:true,json:async()=>JSON.parse(b),arrayBuffer:async()=>b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength)};};
test('geranial scenario retains all exact identities, certificates and bounds',()=>{
  const index=read(folder+'bundle.json'),old=read('docs/data/c17-net-view/bundle.json');
  assert.equal(index.targets.length,6220);
  assert.equal(Object.keys(index.certificate_files).length,2737);
  assert.equal(index.summary.target_status_counts['exact-net-conversion-hypothesis'],2740);
  assert.deepEqual(index.model_summary,report.summary);
  assert.deepEqual(index.forbidden_step_ids,report.forbidden_step_ids);
  assert.deepEqual(index.external_exchange_compound_ids,report.external_exchange_compound_ids);
  assert.deepEqual(index.targets.map(t=>[t.cannabisdb_id,t.compound_id,t.net_status]),report.targets.map(t=>[t.cannabisdb_id,t.compound_id,t.net_status]));
  assert.match(index.medium_annotation_scope,/not recomputed/);
  for(const [cid,ref] of Object.entries(index.certificate_files)){
    const bytes=fs.readFileSync(folder+ref.file);
    assert.equal(bytes.length,ref.bytes);
    assert.equal(createHash('sha256').update(bytes).digest('hex'),ref.sha256);
    if(old.certificate_files[cid]) assert.equal(ref.sha256,old.certificate_files[cid].sha256);
    else assert.deepEqual(JSON.parse(bytes),audit.certificates.find(c=>c.compound_id===cid).certificate);
  }
});
test('all four gain graphs preserve full participants and differentiated evidence',async()=>{
  const ctx={crypto:webcrypto,TextDecoder};vm.runInNewContext(fs.readFileSync('docs/net.js','utf8'),ctx);
  const base=await ctx.NetView.createLoader(fetcher,'geranial-net-view')();
  for(const a of audit.certificates){
    const selected=await base.loadTarget(a.cannabisdb_id), graph=ctx.NetView.project(selected,a.cannabisdb_id);
    assert.equal(graph.steps.length,a.certificate.steps.length);
    for(const s of graph.steps){
      const expected=a.steps.find(x=>x.step_id===s.step_id),edges=graph.edges.filter(e=>e.data.step_id===s.step_id);
      assert.equal(edges.length,s.required_inputs.length*s.outputs.length);
      for(const e of edges){
        assert.deepEqual(Array.from(e.data.certificate_direction_annotation.review_flags),expected.review_flags);
        assert.equal(e.data.certificate_direction_annotation.evidence_class,expected.evidence_class);
        assert.ok(e.data.required_inputs.length&&e.data.outputs.length);
        if(expected.biochemical_evidence) assert.deepEqual(JSON.parse(JSON.stringify(e.data.certificate_direction_annotation.biochemical_evidence)),expected.biochemical_evidence);
      }
      if(s.reaction_id.startsWith('geranial-reduction-hypothesis:')) assert.equal(s.direction_mode,'hypothetical-left-to-right');
    }
  }
});
test('actual viewer mounts geranial target with arrows and readable assay warning',async()=>{
  class Field {constructor(){this.value='';this.textContent='';this.children=[];}addEventListener(e,f){this[e]=f;}replaceChildren(...c){this.children=c;this.value=c[0]?.value||'';}appendChild(c){this.children.push(c);}set innerHTML(v){throw Error('Unsafe HTML');}}
  const fields={};const field=id=>fields[id]||(fields[id]=new Field());
  field('netScope').value='certificates';field('poolHighlight').value='all';
  const cy={items:[],elements(){return {remove:()=>{this.items=[];},removeClass(){}};},nodes(){return {filter(){return {addClass(){}};}};},edges(){return {filter(){return {addClass(){}};}};},add(x){this.items.push(...x);},layout(){return {run(){}};},fit(){},on(){}};
  let style;const fetched=[];
  const ui={crypto:webcrypto,TextDecoder,URLSearchParams,location:{search:'?scenario=geranial&target=CDB000585'},document:{getElementById:field,createElement:()=>new Field()},cytoscape:o=>{style=o.style;return cy;},fetch:async u=>{fetched.push(u);return fetcher(u);}};
  vm.runInNewContext(fs.readFileSync('docs/net.js','utf8'),ui);const app=ui.NetView.mount();await app.load();
  assert.ok(fetched.every(u=>u.startsWith('data/geranial-net-view/')));
  assert.equal(field('netTarget').value,'CDB000585');assert.match(field('netMetrics').textContent,/2740 \/ 6220/);
  assert.equal(style.find(s=>s.selector==='edge').style['target-arrow-shape'],'triangle');
  const edge=cy.items.find(e=>e.data.certificate_direction_annotation?.evidence_class==='substrate-specific-cascade-biochemistry-outside-Cannabis');assert.ok(edge);
  field('netReaction').value=edge.data.step_id;field('netReaction').change();
  assert.match(JSON.stringify(field('netEquation').children),/Substrate-specific cascade evidence outside Cannabis/);
  assert.match(JSON.stringify(field('netEquation').children),/another organism, not Cannabis/);
});
