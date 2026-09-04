const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const script = fs.readFileSync(path.join(__dirname, '../docs/hypotheses.js'), 'utf8');
const context = {module: {exports: {}}};
vm.runInNewContext(script, context);
const {project, createLoader, matchingTargets} = context.module.exports;

test('projects one complete reaction with arrows, coefficients and no omitted cofactors', () => {
  const c = id => ({id, labels: [id], formula: 'X'});
  const h = {id: 'h', reaction_id: 'r', compound_id: 'p', required_inputs: [{compound_id: 'a', coefficient: 2}, {compound_id: 'water', coefficient: 1}], outputs: [{compound_id: 'p', coefficient: 1}, {compound_id: 'water', coefficient: 1}], has_candidate_enzyme_evidence: false};
  const result = project({compounds: ['a', 'water', 'p'].map(c), reaction: {sources: [{source_reaction_id: 'RHEA:1'}]}, hypotheses: [h]}, 'h');
  assert.equal(result.nodes.length, 3);
  assert.equal(result.edges.length, 4);
  assert.equal(result.nodes.find(n => n.data.id === 'water').data.role, 'both');
  assert.equal(result.nodes.find(n => n.data.id === 'a').data.input_coefficient, 2);
  assert.ok(result.nodes.find(n => n.data.id === 'p').data.is_target);
  assert.ok(result.edges.every(e => e.data.required_inputs === h.required_inputs && e.data.outputs === h.outputs));
  assert.ok(result.edges.some(e => e.data.source === 'water' && e.data.target === 'p'));
  assert.ok(result.edges.every(e => e.data.direction_status.includes('hypothetical')));
});

test('loader fetches only requested shards, caches, retries failures and rejects path injection', async () => {
  const rid = 'balanced-equation:' + 'ab'.repeat(32), calls = [];
  let fail = true;
  const loader = createLoader(async url => {calls.push(url); if (fail) {fail = false; return {ok: false, status: 503};} return {ok: true, json: async () => ({[rid]: {reaction: {id: rid}}})};});
  await assert.rejects(loader.reaction(rid), /503/);
  await loader.reaction(rid); await loader.reaction(rid);
  assert.equal(calls.length, 2);
  assert.equal(calls[0], 'data/hypothesis-view/reactions/ab.json');
  assert.throws(() => loader.target('../secret'), /Invalid/);
  await assert.rejects(loader.reaction('../secret'), /Invalid/);
  assert.equal(matchingTargets([{cannabisdb_id: 'CDB000055', label: 'Eugenol'}], 'eUGEnol').length, 1);
});

test('loader revalidates the manifest and versions data requests by published checksum', async () => {
  const calls = [];
  const loader = createLoader(async (url, options) => {calls.push({url, options}); return {ok: true, json: async () => ({files: {'targets/CDB000055.json': {sha256: '1234567890abcdef' + '0'.repeat(48)}}})};});
  await loader.index(); await loader.target('CDB000055');
  assert.equal(calls[0].options.cache, 'no-cache');
  assert.equal(calls[1].url, 'data/hypothesis-view/targets/CDB000055.json?v=1234567890abcdef');
});

test('published bundles render all required participants for every hypothesis', () => {
  const folder = path.join(__dirname, '../docs/data/hypothesis-view');
  const index = JSON.parse(fs.readFileSync(path.join(folder, 'index.json')));
  let count = 0;
  for (const name of Object.keys(index.files).filter(n => n.startsWith('reactions/'))) {
    for (const bundle of Object.values(JSON.parse(fs.readFileSync(path.join(folder, name))))) {
      for (const h of bundle.hypotheses) {
        const graph = project(bundle, h.id);
        const expected = new Set([...h.required_inputs, ...h.outputs].map(m => m.compound_id));
        assert.deepEqual(new Set(graph.nodes.map(n => n.data.id)), expected);
        assert.equal(graph.edges.length, h.required_inputs.length * h.outputs.length);
        assert.equal(new Set(graph.edges.map(e => e.data.reaction_id)).size, 1);
        const screenedIds = (bundle.enzyme_evidence || []).filter(e => (h.evidence_ids || []).includes(e.id)).flatMap(e => (e.screened_proteins || []).map(p => p.accession));
        assert.ok(graph.edges.every(e => screenedIds.every(id => e.data.candidate_protein_ids.includes(id))));
        count++;
      }
    }
  }
  assert.equal(count, index.summary.one_step_hypotheses);
});

test('controls handle gaps, evidence filters and stale responses without displaying the wrong target', async () => {
  class Field {
    constructor() {this.value = ''; this.textContent = ''; this.children = []; this.hidden = false;}
    addEventListener(event, callback) {this[event] = callback;}
    replaceChildren(...children) {this.children = children; this.value = children[0]?.value || '';}
    appendChild(child) {this.children.push(child);}
    set innerHTML(value) {throw new Error('Untrusted content must use textContent');}
  }
  const ids = ['hypothesisCy', 'fit', 'retry', 'message', 'selectedDetails', 'hypothesisSelect', 'targetSelect', 'targetSearch', 'enzymeFilter', 'graphCount', 'equation', 'blockers', 'tests', 'sources', 'evidence', 'hypothesisCount', 'selectionTitle', 'targetState', 'searchCount', 'metrics'];
  const fields = Object.fromEntries(ids.map(id => [id, new Field()])); fields.enzymeFilter.value = 'all';
  const rid = 'balanced-equation:' + 'ab'.repeat(32);
  const targets = [{cannabisdb_id: 'CDB000001', label: 'Eugenol', hypothesis_count: 1, status: 'net-production-hypotheses-found', structure_status: 'known', next_step: 'test'}, {cannabisdb_id: 'CDB000002', label: '<img onerror=alert(1)>', hypothesis_count: 0, status: 'no-exact-encoded-reaction-match', structure_status: 'known', next_step: 'curate'}];
  const h = {id: 'h', reaction_id: rid, compound_id: 'product', required_inputs: [{compound_id: 'input', coefficient: 1}], outputs: [{compound_id: 'product', coefficient: 1}], has_candidate_enzyme_evidence: false, blockers: ['enzyme-unconfirmed'], proposed_tests: ['test activity'], net_target_coefficient: 1};
  const bundle = {reaction: {id: rid, sources: [{source_reaction_id: 'RHEA:1', source_urls: ['javascript:alert(1)', 'https://www.rhea-db.org/rhea/1']}]}, hypotheses: [h], compounds: ['input', 'product'].map(id => ({id, labels: [id]})), enzyme_evidence: []};
  let resolveTarget, resolveReaction;
  const targetPromise = new Promise(resolve => {resolveTarget = resolve;});
  const reactionPromise = new Promise(resolve => {resolveReaction = resolve;});
  const cy = {items: [], elements() {return {remove: () => {this.items = [];}};}, add(items) {this.items.push(...items);}, layout() {return {run() {}};}, fit() {}, on() {}};
  const ui = {module: {exports: {}}, URLSearchParams, location: {search: ''},
    document: {getElementById: id => fields[id], createElement: () => new Field()}, cytoscape: () => cy,
    fetch: async url => ({ok: true, json: async () => url.endsWith('index.json') ? {targets, summary: {carbon_bearing_target_status_counts: {'net-production-hypotheses-found': 1}, carbon_bearing_target_records: 2, cannabisdb_records: 2}} : url.includes('/targets/') ? targetPromise : reactionPromise})};
  vm.runInNewContext(script, ui); ui.module.exports.mount();
  await new Promise(setImmediate);
  fields.enzymeFilter.change(); // Must not cancel target loading.
  resolveTarget({target: targets[0], hypotheses: [{...h, source_reaction_ids: ['RHEA:1']}]});
  await new Promise(setImmediate);
  fields.targetSearch.value = 'CDB000002'; fields.targetSearch.input();
  resolveReaction({[rid]: bundle}); await new Promise(setImmediate);
  assert.equal(cy.items.length, 0);
  assert.match(fields.message.textContent, /No balanced net-production/);
  assert.match(fields.selectionTitle.textContent, /<img/); // Literal text, never markup.
  fields.targetSearch.value = 'Eugenol'; fields.targetSearch.input(); await new Promise(setImmediate);
  assert.equal(cy.items.length, 3);
  assert.match(fields.graphCount.textContent, /1 reaction/);
  assert.equal(fields.sources.children.filter(c => c.children.some(a => a.href?.startsWith('javascript:'))).length, 0);
  fields.enzymeFilter.value = 'candidate'; fields.enzymeFilter.change(); await new Promise(setImmediate);
  assert.equal(cy.items.length, 0);
  assert.match(fields.message.textContent, /No hypotheses match/);
  fields.enzymeFilter.value = 'missing'; fields.enzymeFilter.change(); await new Promise(setImmediate);
  assert.equal(cy.items.length, 3);
});
