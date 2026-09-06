const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const read = name => JSON.parse(fs.readFileSync(path.join(__dirname, '..', name)));
const script = fs.readFileSync(path.join(__dirname, '../docs/net.js'), 'utf8');
const context = {}; vm.runInNewContext(script, context);
const {createLoader, project} = context.NetView;
const bundle = read('docs/data/pg-named-net-view/bundle.json');
const manifest = read('docs/data/pg-named-net-view/index.json');
const net = read('data/reports/phase1-pg-named-net.json');
const fetcher = data => async url => ({ok:true, json:async()=>url.endsWith('index.json')?manifest:data});

test('all three PG comparisons retain exact identities and full directed equations', async()=>{
  const before = JSON.stringify(bundle);
  for (const option of bundle.scenario_options) {
    const loaded = await createLoader(fetcher(bundle), 'pg-named-net-view', option.id)();
    assert.equal(loaded.summary.target_records, 23);
    assert.ok(loaded.view_boundary.startsWith(option.title));
    for (const pair of net.paired_probes) {
      const graph = project(loaded, pair.cannabisdb_id);
      assert.equal(graph.target.compound_id, pair[option.id].compound_id);
      if (option.id !== 'alternative_extended_result') {
        assert.equal(graph.certificate, null);
        assert.equal(graph.nodes.length, 0);
        assert.equal(graph.edges.length, 0);
        continue;
      }
      assert.deepEqual(graph.certificate, pair[option.id]);
      assert.equal(graph.edges.length, graph.steps.reduce((n,s)=>n+s.required_inputs.length*s.outputs.length,0));
      for (const step of graph.steps) {
        assert.ok(!net.forbidden_step_ids.includes(step.step_id));
        const r = net.certificate_reactions.find(r=>r.id===step.reaction_id);
        const forward = step.direction_mode==='hypothetical-left-to-right';
        assert.deepEqual(step.required_inputs, r[forward?'left':'right']);
        assert.deepEqual(step.outputs, r[forward?'right':'left']);
        for (const edge of graph.edges.filter(e=>e.data.step_id===step.step_id)) {
          assert.ok(step.required_inputs.some(p=>p.compound_id===edge.data.source));
          assert.ok(step.outputs.some(p=>p.compound_id===edge.data.target));
          assert.deepEqual(edge.data.required_inputs, step.required_inputs);
          assert.deepEqual(edge.data.outputs, step.outputs);
          if (r.hypothesis_type) {
            assert.equal(edge.data.is_route_sensitivity, true);
            assert.ok(edge.data.hypothesis_assumptions.includes(r.claim_boundary));
          }
        }
      }
    }
  }
  assert.equal(JSON.stringify(bundle), before);
  await assert.rejects(createLoader(fetcher(bundle),'pg-named-net-view','invented')(), /Invalid identity/);
  for (const bad of [{...bundle,scenario_options:[]}, {...bundle,view_scenario:'reaction-first-chemistry'}]) {
    await assert.rejects(createLoader(fetcher(bad),'pg-named-net-view','original_result')(), /Invalid paired/);
  }
});

for (const comparison of ['', 'original_result', 'alternative_baseline_result']) {
  test(`PG controls show correct denominator and empty states: ${comparison || 'default alternative'}`, async()=>{
    class Field {
      constructor(){this.value='';this.textContent='';this.children=[];this.hidden=false;}
      addEventListener(event,callback){this[event]=callback;}
      replaceChildren(...children){this.children=children;this.value=children[0]?.value||'';}
      appendChild(child){this.children.push(child);}
      set innerHTML(value){throw Error('Use textContent');}
    }
    const fields={};
    const field = id => fields[id] ||= new Field();
    field('netScope').value='certificates'; field('poolHighlight').value='all';
    const cy={items:[],elements(){return {remove:()=>{this.items=[];},removeClass(){}};},
      nodes(){return {filter(){return {addClass(){}};}};},edges(){return {filter(){return {addClass(){}};}};},
      add(items){this.items.push(...items);},layout(){return {run(){}};},fit(){},on(){}};
    let options;
    const ui={URLSearchParams,location:{search:'?scenario=pg-named'+(comparison?'&comparison='+comparison:'')},
      document:{getElementById:field,createElement:()=>new Field()},
      cytoscape:value=>{options=value;return cy;},fetch:fetcher(bundle)};
    vm.runInNewContext(script,ui); ui.NetView.mount(); await new Promise(setImmediate);
    const positive=!comparison;
    assert.ok(field('netMetrics').textContent.startsWith((positive?'23':'0')+' / 23'));
    assert.ok(field('netMetrics').textContent.includes('no historical coverage gain'));
    assert.equal(field('netRetry').hidden,true);
    assert.ok(field('netTitle').textContent.includes(positive?'Name-derived':comparison==='original_result'?'Original encoded':'Name-derived'));
    assert.equal(cy.items.length>0,positive);
    assert.ok(options.style.some(s=>s.style?.['target-arrow-shape']==='triangle'));
    if (!positive) {
      assert.equal(field('netScope').value,'all');
      assert.equal(field('netTarget').children.length,23);
    }
  });
}
