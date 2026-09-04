const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const sandbox = {module: {exports: {}}};
vm.runInNewContext(fs.readFileSync(path.join(__dirname, '../docs/phase1-evidence.js'), 'utf8'), sandbox);
const {makeIndex, select} = sandbox.module.exports;

test('joins exact SMARTS variants and counts variants rather than projected edges', () => {
  const a = {reaction_id: 'RHEA:1', reaction_smarts: 'a>>b', evidence_status: 'screened-homology'};
  const b = {...a, reaction_smarts: 'c>>d', evidence_status: 'missing-reference'};
  const index = makeIndex([a, b]);
  const edges = [
    {...a, id: '1', source: 'a', target: 'b'},
    {...a, id: '2', source: 'a', target: 'c'},
    {...b, id: '3', source: 'c', target: 'd'},
    {...a, reaction_id: 'MARTS:1', id: '4', source: 'e', target: 'f'},
    {...a, reaction_smarts: undefined, id: '5', source: 'g', target: 'h'},
  ];
  const selection = select(edges, index, 'screened-homology');
  assert.deepEqual([...selection.edgeIds], ['1', '2']);
  assert.deepEqual([...selection.nodeIds], ['a', 'b', 'c']);
  assert.equal(selection.variants.size, 1);
  assert.equal(select(edges, index, 'all').variants.size, 2);
  assert.equal(select(edges, index, 'no-hits').edgeIds.size, 0);
  assert.throws(() => makeIndex([a, a]), /Duplicate/);
});

test('published overlay covers every expansion variant without inflating edge counts', () => {
  const read = relative => JSON.parse(fs.readFileSync(path.join(__dirname, relative), 'utf8'));
  const report = read('../docs/data/phase1-map-evidence.json');
  const expansion = read('../docs/data/terpene-identity-set-candidate-expansion.json');
  const selection = select(expansion.rows.map((r, i) => ({...r, id: String(i), source: r.precursor_terpene_id, target: r.product_terpene_id})), makeIndex(report.rows), 'all');
  assert.equal(selection.variants.size, report.summary.reaction_variants);
  assert.equal(selection.edgeIds.size, expansion.rows.length);
  assert.deepEqual(read('../data/reports/phase1-map-evidence.json'), report);
});

test('UI adapter updates loaded edges, composes dimming, and shows evidence safely', async () => {
  const row = {reaction_id: 'RHEA:1', reaction_smarts: 'a>>b', evidence_status: 'screened-homology',
    balance_status: 'balanced', search_status: 'hits-found', screened_proteins: ['P1'], core_enzyme_ids: []};
  const fields = {};
  for (const id of ['phase1Status', 'phase1Metrics', 'phase1Matches', 'details']) {
    fields[id] = {value: 'all', disabled: true, textContent: '', children: [],
      addEventListener(name, callback) { this[name] = callback; },
      appendChild(child) { this.children.push(child); }};
  }
  const element = data => ({data: () => data, id: () => data.id, classes: new Set(['faded']),
    toggleClass(name, value) { value ? this.classes.add(name) : this.classes.delete(name); }});
  const edges = [element({...row, id: 'e1', source: 'n1', target: 'n2'}), element({id: 'e2', source: 'n3', target: 'n4'})];
  const nodes = ['n1', 'n2', 'n3', 'n4'].map(id => element({id}));
  const handlers = {};
  const style = {selector() { return this; }, style() { return this; }, update() {}};
  const cy = {edges: () => edges, nodes: () => nodes, batch: callback => callback(), style: () => style,
    on(name, selector, callback) { handlers[name] = callback || selector; }};
  const context = {module: {exports: {}}, queueMicrotask, console,
    document: {getElementById: id => fields[id], createElement: () => ({textContent: ''})},
    fetch: async () => ({ok: true, json: async () => ({rows: [row], summary: {
      reaction_variants: 1, balanced_variants: 1, balance_status_counts: {balanced: 1},
      balanced_variants_with_candidate_enzyme_evidence: 1, balanced_variants_without_candidate_enzyme_evidence: 0}})})};
  vm.runInNewContext(fs.readFileSync(path.join(__dirname, '../docs/phase1-evidence.js'), 'utf8'), context);
  context.module.exports.mount(cy);
  await new Promise(setImmediate);
  assert.equal(fields.phase1Status.disabled, false);
  fields.phase1Status.value = 'screened-homology'; fields.phase1Status.change();
  assert.equal(edges[0].classes.has('phase1-dim'), false);
  assert.equal(edges[1].classes.has('phase1-dim'), true);
  assert.equal(nodes[2].classes.has('phase1-dim'), true);
  assert.equal(edges[0].classes.has('faded'), true); // Existing filters are preserved.
  edges.push(element({...row, id: 'e3', source: 'n3', target: 'n4'}));
  handlers.add(); await new Promise(setImmediate);
  assert.equal(nodes[2].classes.has('phase1-dim'), false);
  handlers.tap({target: edges[0]});
  assert.match(fields.details.children[0].textContent, /not confirmed Cannabis/);
  fields.phase1Status.value = 'all'; fields.phase1Status.change();
  assert.equal(edges[1].classes.has('phase1-dim'), false);
});
