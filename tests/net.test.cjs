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

for (const bundle of [JSON.parse(fs.readFileSync(path.join(__dirname, '../docs/data/net-view/bundle.json'))), sensitivityBundle]) {
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

for (const [bundle, scenario] of [[JSON.parse(fs.readFileSync(path.join(__dirname, '../docs/data/net-view/bundle.json'))), ''], [sensitivityBundle, '?scenario=completions&target=CDB006149']]) {
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
  const cy={items:[],elements(){return {remove:()=>{this.items=[];},removeClass(){}};},nodes(){return {filter(){return {addClass(){}};}};},add(items){this.items.push(...items);},layout(){return {run(){}};},fit(){},on(){}};
  let fail=true;
  const ui={URLSearchParams,location:{search:scenario},document:{getElementById:id=>fields[id],createElement:()=>new Field()},cytoscape:()=>cy,
    fetch:async url=>{if(fail){fail=false;return {ok:false,status:503};}return {ok:true,json:async()=>url.endsWith('index.json')?{file:'bundle.json',sha256:'a'.repeat(64)}:bundle};}};
  vm.runInNewContext(script,ui); const app=ui.NetView.mount(); await new Promise(setImmediate);
  assert.equal(fields.netRetry.hidden,false);
  await app.load(); assert.equal(fields.netRetry.hidden,true); assert.ok(cy.items.length>0);
  if (scenario) {
    assert.equal(fields.netTarget.value,'CDB006149');
    assert.match(fields.netBoundary.textContent,/Completion sensitivity/);
    assert.match(fields.netMessage.textContent,/unverified completion chemistry/);
    assert.ok(cy.items.some(e=>e.data.is_completion_sensitivity));
  }
  const fullCount=cy.items.length;
  fields.poolHighlight.value='pools';fields.poolHighlight.change();assert.equal(cy.items.length,fullCount);
  fields.netReaction.value=fields.netReaction.children[1].value;fields.netReaction.change();
  assert.match(fields.netCounts.textContent,/1 directed reactions shown/);
  assert.ok(fields.netBalance.children.length);assert.match(fields.netEvidence.textContent,/"id"/);
  fields.netScope.value='all';fields.netSearch.value=bundle.targets.find(t=>!t.certificate_compound_id).cannabisdb_id;fields.netSearch.input();
  assert.equal(cy.items.length,0);assert.match(fields.netMessage.textContent,/No net-conversion certificate/);
  fields.netSearch.value='no such compound xyz';fields.netSearch.input();assert.equal(cy.items.length,0);assert.equal(fields.netTarget.disabled,true);
});
}
