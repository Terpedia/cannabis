const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const script = fs.readFileSync(path.join(__dirname, '../docs/net.js'), 'utf8');
const context = {}; vm.runInNewContext(script, context);
const {project, matchingTargets, createLoader} = context.NetView;
const bundle = JSON.parse(fs.readFileSync(path.join(__dirname, '../docs/data/net-view/bundle.json')));
const sensitivityBundle = JSON.parse(fs.readFileSync(path.join(__dirname, '../docs/data/completion-net-view/bundle.json')));
const catalogBundle = JSON.parse(fs.readFileSync(path.join(__dirname, '../docs/data/catalog-net-view/bundle.json')));
const catalogEvidence = JSON.parse(fs.readFileSync(path.join(__dirname, '../docs/data/catalog-net-view/evidence.json')));
const updatedCatalogBundle = context.NetView.applyEvidence(catalogBundle,catalogEvidence);

test('catalog supplement loader preserves chemistry, fails closed, and retains original gaps', async()=>{
  const manifest=JSON.parse(fs.readFileSync(path.join(__dirname,'../docs/data/catalog-net-view/index.json')));
  const before=JSON.stringify(catalogBundle),calls=[];
  const fetcher=async url=>{calls.push(url);return {ok:true,json:async()=>url.endsWith('index.json')?manifest:url.includes('evidence.json')?catalogEvidence:catalogBundle};};
  const updated=await createLoader(fetcher,'catalog-net-view')();
  assert.equal(calls.length,3);
  assert.equal(JSON.stringify(catalogBundle),before);
  assert.equal(updated.certificates.length,catalogBundle.certificates.length);
  for(const [i,c] of updated.certificates.entries()) {
    for(const k of ['steps','net_exports','external_net_consumption','zero_net_internal_participants']) assert.equal(c[k],catalogBundle.certificates[i][k]);
  }
  assert.equal(updated.reactions.filter(r=>r.is_new_catalog_candidate).length,97);
  assert.equal(updated.reactions.filter(r=>r.missing_candidate_evidence).length,368);
  assert.equal(matchingTargets(updated.targets,'','enzyme-gaps').length,202);
  const uric=project(updated,'CDB004839');
  assert.equal(uric.certificate.missing_candidate_reaction_ids.length,0);
  assert.ok(uric.certificate.baseline_missing_candidate_reaction_ids.length);
  assert.ok(uric.steps.every(s=>s.evidence.length));
  assert.ok(uric.edges.some(e=>e.data.is_new_catalog_candidate));
  const bad=JSON.parse(JSON.stringify(catalogEvidence));bad.certificate_updates[0].missing_candidate_reaction_ids=[];
  assert.throws(()=>context.NetView.applyEvidence(catalogBundle,bad),/mismatch/);
  await assert.rejects(createLoader(async url=>url.includes('evidence.json')?{ok:false,status:503}:fetcher(url),'catalog-net-view')(),/503/);
  const wrong=JSON.parse(JSON.stringify(catalogEvidence));wrong.source_sha256['data/reports/phase1-catalog-net-gaps.json']='0'.repeat(64);
  await assert.rejects(createLoader(async url=>url.includes('evidence.json')?{ok:true,json:async()=>wrong}:fetcher(url),'catalog-net-view')(),/snapshot mismatch/);
});

for (const bundle of [JSON.parse(fs.readFileSync(path.join(__dirname, '../docs/data/net-view/bundle.json'))), sensitivityBundle, catalogBundle, updatedCatalogBundle]) {
test(`every ${bundle.view_scenario || 'baseline'} net certificate projects all participants, coefficients, extents and candidate evidence`, () => {
  let count = 0;
  for (const target of bundle.targets.filter(t => t.certificate_compound_id)) {
    const graph = project(bundle, target.cannabisdb_id);
    assert.equal(graph.steps.length, graph.certificate.steps.length);
    const participants = new Set(); let arrows = 0;
    for (const step of graph.steps) {
      const source = bundle.reactions.find(r => r.id === step.reaction_id);
      const forward = step.direction_mode === 'hypothetical-left-to-right';
      assert.equal(step.required_inputs, source[forward ? 'left' : 'right']);
      assert.equal(step.outputs, source[forward ? 'right' : 'left']);
      [...step.required_inputs, ...step.outputs].forEach(m => participants.add(m.compound_id));
      arrows += step.required_inputs.length * step.outputs.length;
      const edges = graph.edges.filter(e => e.data.step_id === step.step_id);
      const proteins = step.evidence.flatMap(e => (e.screened_proteins || []).map(p => p.accession));
      for (const edge of edges) {
        assert.equal(edge.data.extent, step.extent);
        assert.equal(edge.data.required_inputs, step.required_inputs);
        assert.equal(edge.data.outputs, step.outputs);
        assert.equal(edge.data.is_completion_sensitivity, !!source.is_completion_sensitivity);
        assert.equal(edge.data.missing_candidate_evidence, !!source.missing_candidate_evidence);
        if(source.missing_candidate_evidence) assert.equal(step.evidence.length, 0);
        assert.ok(proteins.every(p => edge.data.candidate_protein_ids.includes(p)));
      }
      const selected = project(bundle, target.cannabisdb_id, step.step_id);
      assert.equal(selected.steps.length, 1);
      assert.equal(selected.edges.length, edges.length);
      assert.equal(selected.certificate, graph.certificate);
    }
    assert.deepEqual(new Set(graph.nodes.map(n => n.data.id)), participants);
    assert.equal(graph.edges.length, arrows);
    assert.equal(new Set(graph.edges.map(e=>e.data.id)).size, arrows);
    for (const node of graph.nodes) assert.equal(node.data.is_pool, graph.certificate.zero_net_internal_participants.includes(node.data.id));
    count++;
  }
  assert.equal(count, bundle.summary.target_status_counts['exact-net-conversion-hypothesis']);
  const gap = bundle.targets.find(t=>!t.certificate_compound_id);
  assert.equal(project(bundle, gap.cannabisdb_id).nodes.length, 0);
  assert.equal(matchingTargets(bundle.targets, '', 'all').length, bundle.targets.length);
  assert.equal(matchingTargets(bundle.targets, '', 'certificates').length, count);
  assert.equal(matchingTargets(bundle.targets, '', 'enzyme-gaps').length, bundle.targets.filter(t=>t.certificate_compound_id && t.missing_candidate_reaction_ids?.length).length);
  assert.throws(()=>project(bundle,'missing'), /Unknown target/);
});
}

test('loader revalidates manifest, versions bundles, rejects external paths and allows retry', async () => {
  const calls = []; let failed = false;
  const loader = createLoader(async (url, options) => {
    calls.push({url, options});
    if (!failed) {failed = true; return {ok:false,status:503};}
    return {ok:true,json:async()=>url.endsWith('index.json') ? {file:'bundle.json',sha256:'a'.repeat(64)} : bundle};
  });
  await assert.rejects(loader(), /503/); assert.equal(await loader(), bundle);
  assert.equal(calls[1].options.cache,'no-cache');
  assert.equal(calls[2].url,'data/net-view/bundle.json?v='+'a'.repeat(16));
  await assert.rejects(createLoader(async()=>({ok:true,json:async()=>({file:'https://external.invalid/file',sha256:'a'.repeat(64)})}))(), /Invalid/);
  assert.throws(()=>createLoader(()=>{}, '../external'), /Invalid scenario/);
  const sensitivityCalls=[];
  await createLoader(async url=>{sensitivityCalls.push(url);return {ok:true,json:async()=>url.endsWith('index.json')?{file:'bundle.json',sha256:'b'.repeat(64)}:sensitivityBundle};},'completion-net-view')();
  assert.deepEqual(sensitivityCalls,['data/completion-net-view/index.json','data/completion-net-view/bundle.json?v='+'b'.repeat(16)]);
});

for (const [bundle, scenario] of [[JSON.parse(fs.readFileSync(path.join(__dirname, '../docs/data/net-view/bundle.json'))), ''], [sensitivityBundle, '?scenario=completions&target=CDB006149'], [catalogBundle, '?scenario=catalog&target=CDB006137'], [updatedCatalogBundle, '?scenario=catalog&target=CDB006137']]) {
test(`controls ${scenario || 'baseline'} retain full balances, clear gaps, highlight without hiding, and recover from load errors`, async () => {
  class Field {
    constructor(){this.value='';this.textContent='';this.children=[];this.hidden=false;}
    addEventListener(event, callback){this[event]=callback;}
    replaceChildren(...children){this.children=children;this.value=children[0]?.value||'';}
    appendChild(child){this.children.push(child);}
    set innerHTML(value){throw Error('Use safe textContent');}
  }
  const ids=['netCy','netFit','poolHighlight','netSearch','netScope','netTarget','netReaction','netRetry','netCounts','netBalance','netEquation','netSources','netEvidence','netDetails','netTitle','netStatus','netMessage','netMetrics','netMatches','netBoundary'];
  const fields=Object.fromEntries(ids.map(id=>[id,new Field()]));fields.netScope.value='certificates';fields.poolHighlight.value='all';
  const cy={items:[],muted:[],elements(){return {remove:()=>{this.items=[];},removeClass:()=>{this.muted=[];}};},nodes(){return {filter(){return {addClass(){}};}};},edges(){return {filter:predicate=>({addClass:()=>{this.muted=this.items.filter(e=>e.data.source && predicate({data:key=>e.data[key]}));}})};},add(items){this.items.push(...items);},layout(){return {run(){}};},fit(){},on(){}};
  let fail=true, cyOptions; const fetched=[];
  const ui={URLSearchParams,location:{search:scenario},document:{getElementById:id=>fields[id],createElement:()=>new Field()},cytoscape:options=>{cyOptions=options;return cy;},
    fetch:async url=>{fetched.push(url);if(fail){fail=false;return {ok:false,status:503};}return {ok:true,json:async()=>url.endsWith('index.json')?{file:'bundle.json',sha256:'a'.repeat(64)}:bundle};}};
  vm.runInNewContext(script,ui); const app=ui.NetView.mount(); await new Promise(setImmediate);
  assert.equal(fields.netRetry.hidden,false);
  await app.load(); assert.equal(fields.netRetry.hidden,true); assert.ok(cy.items.length>0);
  if (scenario.includes('completions')) {
    assert.equal(fields.netTarget.value,'CDB006149');
    assert.match(fields.netBoundary.textContent,/Completion sensitivity/);
    assert.match(fields.netMessage.textContent,/unverified completion chemistry/);
    assert.ok(cy.items.some(e=>e.data.is_completion_sensitivity));
  }
  if (scenario.includes('catalog')) {
    assert.equal(fetched.at(-1),'data/catalog-net-view/bundle.json?v='+'a'.repeat(16));
    assert.equal(cyOptions.style.find(s=>s.selector==='edge').style['target-arrow-shape'],'triangle');
    const gapStyle=cyOptions.style.find(s=>s.selector==='edge[?missing_candidate_evidence]').style;
    assert.equal(gapStyle['line-color'],'#ff7777');
    assert.equal(gapStyle['target-arrow-color'],'#ff7777');
    assert.match(fields.netBoundary.textContent,/Chemistry-only|Balanced-catalog/);
    assert.match(fields.netMetrics.textContent,/304 \/ 6220.*enzyme gaps included/);
    assert.ok(cy.items.some(e=>e.data.missing_candidate_evidence));
    const originalCount=cy.items.length;
    fields.poolHighlight.value='enzyme-gaps';fields.poolHighlight.change();
    assert.equal(cy.items.length,originalCount);
    assert.ok(cy.muted.every(e=>!e.data.missing_candidate_evidence));
    assert.equal(cy.muted.length,cy.items.filter(e=>e.data.source && !e.data.missing_candidate_evidence).length);
    fields.netScope.value='enzyme-gaps';fields.netScope.change();
    assert.match(fields.netMatches.textContent,bundle.evidence_summary ? /202 matches/ : /203 matches/);
    const selected=project(bundle,fields.netTarget.value);
    fields.netReaction.value=selected.steps.find(s=>!s.evidence.length).step_id;fields.netReaction.change();
    assert.match(fields.netEvidence.textContent,/No candidate enzyme evidence/);
    fields.netReaction.value='';fields.netReaction.change();
  }
  const fullCount=cy.items.length;
  fields.poolHighlight.value='pools';fields.poolHighlight.change();assert.equal(cy.items.length,fullCount);
  fields.netReaction.value=project(bundle,fields.netTarget.value).steps.find(s=>s.evidence.length).step_id;fields.netReaction.change();
  assert.match(fields.netCounts.textContent,/1 directed reactions shown/);
  assert.ok(fields.netBalance.children.length);assert.match(fields.netEvidence.textContent,/"id"/);
  fields.netScope.value='all';fields.netSearch.value=bundle.targets.find(t=>!t.certificate_compound_id).cannabisdb_id;fields.netSearch.input();
  assert.equal(cy.items.length,0);assert.match(fields.netMessage.textContent,/No net-conversion certificate/);
  fields.netSearch.value='no such compound xyz';fields.netSearch.input();assert.equal(cy.items.length,0);assert.equal(fields.netTarget.disabled,true);
});
}
