const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const {webcrypto,createHash}=require('node:crypto');
const context={crypto:webcrypto,TextDecoder};
vm.runInNewContext(fs.readFileSync('docs/net.js','utf8'),context);
const read=p=>JSON.parse(fs.readFileSync(p));
const report=read('data/reports/phase1-c17-elongation-net.json');
const alkane=read('data/reports/phase1-alkane-net.json');
const audit=read('data/reports/phase1-c17-evidence-audit.json');
const folder='docs/data/c17-net-view/';
test('C17 controls select the correct scenario and display direction warnings without hiding inputs',async()=>{
  class Field {
    constructor(){this.value='';this.textContent='';this.children=[];}
    addEventListener(event,callback){this[event]=callback;}
    replaceChildren(...children){this.children=children;this.value=children[0]?.value||'';}
    appendChild(child){this.children.push(child);}
    set innerHTML(value){throw Error('Unsafe HTML');}
  }
  const ids=['netCy','netFit','poolHighlight','netSearch','netScope','netTarget','netReaction','netRetry','netCounts','netBalance','netEquation','netSources','netEvidence','netDetails','netTitle','netStatus','netMessage','netMetrics','netMatches','netBoundary'];
  const fields=Object.fromEntries(ids.map(id=>[id,new Field()]));
  fields.netScope.value='certificates';fields.poolHighlight.value='all';
  const cy={items:[],muted:[],elements(){return {remove:()=>{this.items=[];},removeClass:()=>{this.muted=[];}};},nodes(){return {filter(){return {addClass(){}};}};},edges(){return {filter:predicate=>({addClass:()=>{this.muted=this.items.filter(e=>e.data.source&&predicate({data:k=>e.data[k]}));}})};},add(items){this.items.push(...items);},layout(){return {run(){}};},fit(){},on(){}};
  const fetched=[];let style;
  const ui={crypto:webcrypto,TextDecoder,URLSearchParams,location:{search:'?scenario=c17&target=CDB000155'},
    document:{getElementById:id=>fields[id],createElement:()=>new Field()},
    cytoscape:opts=>{style=opts.style;return cy;},fetch:async url=>{
      fetched.push(url);const bytes=fs.readFileSync('docs/'+url.split('?')[0]);
      return {ok:true,json:async()=>JSON.parse(bytes),arrayBuffer:async()=>bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength)};
    }};
  vm.runInNewContext(fs.readFileSync('docs/net.js','utf8'),ui);
  const app=ui.NetView.mount();await app.load();
  assert.ok(fetched.every(u=>u.startsWith('data/c17-net-view/')));
  assert.equal(fields.netTarget.value,'CDB000155');
  assert.match(fields.netMetrics.textContent,/2736 \/ 6220/);
  assert.equal(style.find(s=>s.selector==='edge').style['target-arrow-shape'],'triangle');
  const edge=cy.items.find(e=>e.data.certificate_direction_annotation?.review_flags.includes('peroxide-consuming-oxygen-producing-direction-review'));
  assert.ok(edge);
  const count=cy.items.length;fields.poolHighlight.value='direction-review';fields.poolHighlight.change();
  assert.equal(cy.items.length,count);assert.ok(cy.muted.length);
  fields.netReaction.value=edge.data.step_id;fields.netReaction.change();
  assert.match(JSON.stringify(fields.netEquation.children),/consumes peroxide and produces oxygen/);
  assert.match(JSON.stringify(fields.netEquation.children),/does not establish a Cannabis enzyme/);
});
test('C17 map retains all 6220 targets, 2733 hashed certificates and prior witnesses',()=>{
  const index=read(folder+'bundle.json'), old=read('docs/data/ketone-stereo-net-view/bundle.json');
  assert.equal(index.targets.length,6220);
  assert.equal(Object.keys(index.certificate_files).length,2733);
  assert.equal(index.summary.target_status_counts['exact-net-conversion-hypothesis'],2736);
  assert.deepEqual(index.model_summary,report.summary);
  assert.deepEqual(index.forbidden_step_ids,report.forbidden_step_ids);
  assert.deepEqual(index.external_exchange_compound_ids,report.external_exchange_compound_ids);
  assert.deepEqual(index.targets.map(t=>[t.cannabisdb_id,t.compound_id,t.net_status]),report.targets.map(t=>[t.cannabisdb_id,t.compound_id,t.net_status]));
  assert.match(index.medium_annotation_scope,/not recomputed/);
  assert.match(index.view_boundary,/not confirmed Cannabis/);
  const newCerts=[...alkane.new_certificates,...report.new_certificates];
  for(const [cid,ref] of Object.entries(index.certificate_files)){
    const bytes=fs.readFileSync(folder+ref.file);
    assert.equal(bytes.length,ref.bytes);
    assert.equal(createHash('sha256').update(bytes).digest('hex'),ref.sha256);
    if(old.certificate_files[cid]) assert.equal(ref.sha256,old.certificate_files[cid].sha256);
    else assert.deepEqual(JSON.parse(bytes),newCerts.find(c=>c.compound_id===cid));
  }
});
test('six new projected maps retain every input, full certificate and direction-specific warnings',async()=>{
  const fetcher=async url=>{const bytes=fs.readFileSync('docs/'+url.split('?')[0]);return {ok:true,json:async()=>JSON.parse(bytes),arrayBuffer:async()=>bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength)};};
  const base=await context.NetView.createLoader(fetcher,'c17-net-view')();
  for(const cert of [...alkane.new_certificates,...report.new_certificates]){
    const t=report.targets.find(t=>t.compound_id===cert.compound_id);
    const selected=await base.loadTarget(t.cannabisdb_id);
    const graph=context.NetView.project(selected,t.cannabisdb_id);
    assert.equal(graph.steps.length,cert.steps.length);
    const expected=audit.certificates.find(c=>c.cannabisdb_id===t.cannabisdb_id);
    for(const s of graph.steps){
      const edges=graph.edges.filter(e=>e.data.step_id===s.step_id);
      assert.equal(edges.length,s.required_inputs.length*s.outputs.length);
      assert.ok(edges.every(e=>e.data.required_inputs.length&&e.data.outputs.length));
      if(expected){
        const a=expected.steps.find(a=>a.step_id===s.step_id);
        assert.ok(a);
        for(const e of edges){
          assert.deepEqual(Array.from(e.data.certificate_direction_annotation.review_flags),a.review_flags);
          assert.equal(e.data.certificate_direction_annotation.evidence_class,a.evidence_class);
          if(a.review_flags.length) assert.ok(e.data.direction_review_id);
        }
      }
      if(s.reaction_id.startsWith('odd-chain-elongation-hypothesis:')){
        assert.equal(s.direction_mode,'hypothetical-left-to-right');
        assert.equal(s.reaction.is_route_sensitivity,true);
        assert.equal(s.reaction.enzyme_evidence_ids.length,0);
      }
    }
  }
});
