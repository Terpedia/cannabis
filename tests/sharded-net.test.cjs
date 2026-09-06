const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const {webcrypto}=require('node:crypto');
const root=path.join(__dirname,'../docs');
const script=fs.readFileSync(path.join(root,'net.js'),'utf8');
const index=JSON.parse(fs.readFileSync(path.join(root,'data/local-speciation-net-view/bundle.json')));
const model=JSON.parse(fs.readFileSync(path.join(__dirname,'../data/reports/phase1-local-speciation-net.json')));
const context={crypto:webcrypto,TextDecoder};vm.runInNewContext(script,context);
const plain=x=>JSON.parse(JSON.stringify(x));
function response(bytes){return {ok:true,json:async()=>JSON.parse(bytes),arrayBuffer:async()=>bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength)};}
function fetchFiles(calls){return async url=>{calls.push(url);const file=url.split('?')[0];assert.ok(file.startsWith('data/local-speciation-net-view/'));return response(fs.readFileSync(path.join(root,file)));};}

test('on-demand view retains all 6220 targets, fetches no chemistry for gaps, and projects every new certificate',async()=>{
  const calls=[];const load=context.NetView.createLoader(fetchFiles(calls),'local-speciation-net-view');
  const base=await load();assert.equal(calls.length,2);assert.equal(base.targets.length,6220);
  assert.equal(base.summary.target_status_counts['exact-net-conversion-hypothesis'],2723);
  assert.equal(base.certificates.length,0);assert.equal(base.reactions.length,0);
  const gap=base.targets.find(t=>!t.certificate_compound_id);
  const empty=await base.loadTarget(gap.cannabisdb_id);
  assert.equal(calls.length,2);assert.equal(context.NetView.project(empty,gap.cannabisdb_id).certificate,null);
  const expected=new Map(model.new_certificates.map(c=>[c.compound_id,c]));
  for(const target of base.targets.filter(t=>t.new_certificate)){
    const selected=await base.loadTarget(target.cannabisdb_id);
    const graph=context.NetView.project(selected,target.cannabisdb_id);
    assert.deepEqual(plain(graph.certificate),expected.get(target.compound_id));
    assert.equal(selected.certificates.length,1);
    assert.equal(graph.edges.length,graph.steps.reduce((n,s)=>n+s.required_inputs.length*s.outputs.length,0));
    for(const step of graph.steps){
      assert.ok(!base.forbidden_step_ids.includes(step.step_id));
      if(step.reaction.hypothesis_type==='local-speciation'){
        assert.equal(step.reaction.is_route_sensitivity,true);
        assert.ok(step.reaction.hypothesis_assumptions.includes(step.reaction.claim_boundary));
      }
    }
  }
  assert.equal(calls.filter(u=>u.includes('/chemistry.json')).length,1);
  assert.ok(calls.filter(u=>u.includes('/certificates/')).length<=85);
  assert.equal(context.NetView.matchingTargets(base.targets,'','enzyme-gaps').length,base.targets.filter(t=>t.certificate_compound_id&&t.missing_candidate_reaction_count).length);
  await assert.rejects(base.loadTarget('invented'),/Unknown target/);
});

test('corrupt certificate fails hash verification and can be retried without poisoning chemistry cache',async()=>{
  const calls=[];const original=fetchFiles(calls);let corrupt=true;
  const base=await context.NetView.createLoader(async url=>{
    if(url.includes('/certificates/')&&corrupt){corrupt=false;const b=fs.readFileSync(path.join(root,url.split('?')[0]));b[0]^=1;return response(b);}
    return original(url);
  },'local-speciation-net-view')();
  const target=base.targets.find(t=>t.new_certificate);
  await assert.rejects(base.loadTarget(target.cannabisdb_id),/checksum mismatch/);
  const good=await base.loadTarget(target.cannabisdb_id);
  assert.equal(good.certificates[0].compound_id,target.compound_id);
  assert.equal(calls.filter(u=>u.includes('/chemistry.json')).length,1);
});

test('unsafe paths and missing certificate references fail closed before shard fetches',async()=>{
  for(const kind of ['path','identity','count']){
    const bad=structuredClone(index);
    if(kind==='path')bad.shared_chemistry.file='../outside.json';
    if(kind==='identity')delete bad.certificate_files[bad.targets.find(t=>t.certificate_compound_id).compound_id];
    if(kind==='count')bad.targets[0].missing_candidate_reaction_count=-1;
    const calls=[];const original=fetchFiles(calls);
    await assert.rejects(context.NetView.createLoader(async url=>url.includes('/bundle.json')?response(Buffer.from(JSON.stringify(bad))):original(url),'local-speciation-net-view')(),/Invalid/);
    assert.ok(calls.every(u=>!u.includes('/chemistry.json')&&!u.includes('/certificates/')));
  }
});

test('target switches and empty searches cannot render a stale download',async()=>{
  class Field{
    constructor(){this.value='';this.textContent='';this.children=[];this.hidden=false;}
    addEventListener(e,f){this[e]=f;}
    replaceChildren(...c){this.children=c;this.value=c[0]?.value||'';}
    appendChild(c){this.children.push(c);}
    set innerHTML(v){throw Error('Use textContent');}
  }
  const fields={};const field=id=>fields[id] ||= new Field();
  field('netScope').value='all';field('poolHighlight').value='all';
  const cy={items:[],elements(){return{remove:()=>{this.items=[];},removeClass(){}};},
    nodes(){return{filter(){return{addClass(){}};}};},edges(){return{filter(){return{addClass(){}};}};},
    add(x){this.items.push(...x);},layout(){return{run(){}};},fit(){},on(){}};
  const calls=[],original=fetchFiles(calls);let release;let pending=false;
  const first=index.targets.find(t=>t.new_certificate),gap=index.targets.find(t=>!t.certificate_compound_id);
  const ui={crypto:webcrypto,TextDecoder,URLSearchParams,location:{search:'?scenario=speciation&target='+first.cannabisdb_id},
    document:{getElementById:field,createElement:()=>new Field()},cytoscape:()=>cy,
    fetch:async url=>{if(url.includes('/certificates/')&&!pending){pending=true;await new Promise(r=>{release=r;});}return original(url);}};
  vm.runInNewContext(script,ui);ui.NetView.mount();
  for(let i=0;i<50&&!release;i++)await new Promise(setImmediate);
  assert.ok(release);assert.equal(cy.items.length,0);
  field('netTarget').value=gap.cannabisdb_id;await field('netTarget').change();
  assert.ok(field('netTitle').textContent.includes(gap.cannabisdb_id));
  release();await new Promise(r=>setTimeout(r,100));
  assert.equal(cy.items.length,0);assert.ok(field('netTitle').textContent.includes(gap.cannabisdb_id));
  field('netTarget').value=first.cannabisdb_id;await field('netTarget').change();
  assert.ok(cy.items.length>0);assert.ok(field('netTitle').textContent.includes(first.cannabisdb_id));
  field('netSearch').value='no such exact search token';field('netSearch').input();
  assert.equal(cy.items.length,0);assert.equal(field('netTitle').textContent,'No matching target');
  assert.ok(field('netMetrics').textContent.startsWith('2723 / 6220'));
});
