(function(root) {
  'use strict';
  const base = 'data/hypothesis-view/';
  function matchingTargets(targets, query) {
    const text = query.trim().toLowerCase();
    return targets.filter(t => `${t.cannabisdb_id} ${t.label}`.toLowerCase().includes(text));
  }
  function project(bundle, hypothesisId, compact = false) {
    const h = bundle.hypotheses.find(row => row.id === hypothesisId);
    if (!h) throw new Error('Hypothesis missing from reaction shard');
    const compounds = new Map(bundle.compounds.map(c => [c.id, c]));
    const inputs = new Map(h.required_inputs.map(m => [m.compound_id, m.coefficient]));
    const outputs = new Map(h.outputs.map(m => [m.compound_id, m.coefficient]));
    const ids = [...new Set([...inputs.keys(), ...outputs.keys()])];
    const groups = [ids.filter(id => inputs.has(id) && !outputs.has(id)), ids.filter(id => inputs.has(id) && outputs.has(id)), ids.filter(id => !inputs.has(id) && outputs.has(id))];
    let offset = 0;
    const offsets = groups.map(group => { const start = offset; if (group.length) offset += Math.ceil(group.length / 2) * 110 + 70; return start; });
    const nodes = groups.flatMap((group, column) => group.map((id, row) => {
      const c = compounds.get(id);
      if (!c) throw new Error('Required compound is missing');
      const label = c.labels?.[0] || `${c.formula} · ${id.split(':')[1].slice(0, 6)}`;
      return {data: {id, label, kind: 'compound', compound: c, is_target: id === h.compound_id,
        input_coefficient: inputs.get(id) || 0, output_coefficient: outputs.get(id) || 0,
        role: column === 0 ? 'input' : column === 1 ? 'both' : 'output'},
        position: compact ? {x: (row % 2) * 160, y: offsets[column] + Math.floor(row / 2) * 110} : {x: column * 230, y: (row - (group.length - 1) / 2) * 110}};
    }));
    const sourceIds = [...new Set(bundle.reaction.sources.map(s => s.source_reaction_id))];
    const candidateProteinIds = [...new Set((bundle.enzyme_evidence || []).filter(e => (h.evidence_ids || []).includes(e.id)).flatMap(e => [
      ...(e.screened_proteins || []).map(p => p.accession), ...(e.enzyme_ids || []),
      ...(e.enzyme_evidence?.screened_homology_proteins || []), ...(e.enzyme_evidence?.direction_unresolved_family_proteins || [])
    ]))];
    const edges = h.required_inputs.flatMap((input, i) => h.outputs.map((output, j) => ({data: {
      id: `${h.id}:${i}:${j}`, source: input.compound_id, target: output.compound_id,
      label: sourceIds[0] || 'reaction hypothesis', kind: 'reaction-projection',
      hypothesis_id: h.id, reaction_id: h.reaction_id,
      has_candidate_enzyme_evidence: h.has_candidate_enzyme_evidence,
      enzyme_evidence_ids: h.evidence_ids || [],
      candidate_protein_ids: candidateProteinIds,
      direction_status: 'hypothetical; physiological direction unestablished',
      required_inputs: h.required_inputs, outputs: h.outputs,
      claim_boundary: 'Projection of one full reaction; every input is required. Not an independent reaction or an atom-flow edge.'
    }})));
    return {hypothesis: h, nodes, edges};
  }
  function createLoader(fetcher) {
    const cache = new Map();
    let files = {};
    async function load(relative) {
      const version = files[relative]?.sha256;
      const url = base + relative + (version ? '?v=' + version.slice(0, 16) : '');
      if (!cache.has(relative)) cache.set(relative, fetcher(url, {cache: relative === 'index.json' ? 'no-cache' : 'default'}).then(response => {
        if (!response.ok) throw new Error(`Could not load ${relative} (HTTP ${response.status})`);
        return response.json();
      }).catch(error => { cache.delete(relative); throw error; }));
      return cache.get(relative);
    }
    return {index: async () => { const index = await load('index.json'); files = index.files || {}; return index; }, target: id => {
      if (!/^CDB\d+$/.test(id)) throw new Error('Invalid CannabisDB identifier');
      return load(`targets/${id}.json`);
    }, reaction: async id => {
      if (!/^balanced-equation:[a-f0-9]{64}$/.test(id)) throw new Error('Invalid reaction identifier');
      const shard = await load(`reactions/${id.split(':')[1].slice(0, 2)}.json`);
      if (!shard[id]) throw new Error('Reaction missing from shard');
      return shard[id];
    }};
  }
  function mount() {
    const $ = id => document.getElementById(id);
    const loader = createLoader(url => fetch(url));
    let index, targetData, generation = 0, retryAction;
    const cy = cytoscape({container: $('hypothesisCy'), elements: [], layout: {name: 'preset'}, style: [
      {selector: 'node', style: {'label': 'data(label)', 'background-color': '#77b9ef', 'color': '#e8eef8', 'font-size': 16, 'text-wrap': 'wrap', 'text-max-width': 120, 'text-valign': 'bottom', 'text-margin-y': 8, 'width': 48, 'height': 48, 'border-width': 2, 'border-color': '#bfd9ef'}},
      {selector: 'node[role="output"]', style: {'background-color': '#65d6a0'}},
      {selector: 'node[role="both"]', style: {'background-color': '#f2bd65'}},
      {selector: 'node[?is_target]', style: {'border-width': 5, 'border-color': '#ffffff'}},
      {selector: 'edge', style: {'line-color': '#c78cff', 'target-arrow-color': '#c78cff', 'target-arrow-shape': 'triangle', 'curve-style': 'bezier', 'line-style': 'dashed', 'width': 2}},
    ]});
    cy.on('tap', 'node,edge', event => { $('selectedDetails').textContent = JSON.stringify(event.target.data(), null, 2); });
    $('fit').addEventListener('click', () => cy.fit(cy.elements(), 45));
    function options(select, rows) {
      select.replaceChildren(...rows.map(([value, label]) => { const option = document.createElement('option'); option.value = value; option.textContent = label; return option; }));
    }
    function section(id, heading, lines) {
      const title = document.createElement('h2'); title.textContent = heading;
      const list = document.createElement('ul');
      for (const line of lines) { const item = document.createElement('li'); item.textContent = line; list.appendChild(item); }
      $(id).replaceChildren(title, list);
    }
    function clearGraph() {
      cy.elements().remove();
      for (const id of ['equation', 'blockers', 'tests', 'sources', 'evidence', 'graphCount', 'selectedDetails']) $(id).textContent = '';
    }
    function error(error, action, token) {
      if (token !== generation) return;
      clearGraph(); $('message').textContent = error.message; $('retry').hidden = false; retryAction = action;
    }
    $('retry').addEventListener('click', () => { $('retry').hidden = true; retryAction?.(); });
    async function selectHypothesis() {
      const token = ++generation;
      clearGraph(); $('retry').hidden = true;
      const row = targetData?.hypotheses.find(h => h.id === $('hypothesisSelect').value);
      if (!row) { $('message').textContent = 'No hypotheses match this evidence filter.'; return; }
      $('message').textContent = 'Loading selected reaction…';
      try {
        const bundle = await loader.reaction(row.reaction_id);
        if (token !== generation) return;
        const graph = project(bundle, row.id, $('hypothesisCy').clientWidth < 600), h = graph.hypothesis;
        cy.add([...graph.nodes, ...graph.edges]); cy.layout({name: 'preset', fit: true, padding: 30}).run();
        $('graphCount').textContent = `${graph.nodes.length} compounds · 1 reaction · ${graph.edges.length} projected edges`;
        $('message').textContent = `Balanced equation · ${h.has_candidate_enzyme_evidence ? 'candidate enzyme evidence attached' : 'no candidate enzyme evidence attached'} · pathway blocked`;
        const labels = new Map(graph.nodes.map(n => [n.data.id, n.data.label]));
        const format = members => members.map(m => `${m.coefficient} × ${labels.get(m.compound_id)}`).join(' + ');
        section('equation', 'Complete equation — hypothetical direction', [format(h.required_inputs) + ' → ' + format(h.outputs), `Net target produced: ${h.net_target_coefficient}. Every input remains unavailable or unestablished in this pathway assessment.`]);
        section('blockers', 'Unresolved requirements', h.blockers);
        section('tests', 'Proposed tests', h.proposed_tests);
        section('sources', 'Reaction sources', bundle.reaction.sources.map(s => `${s.source_reaction_id} · ${s.source_layer} · source left corresponds to canonical ${s.source_left_corresponds_to}`));
        for (const url of new Set(bundle.reaction.sources.flatMap(s => s.source_urls || []))) {
          if (typeof url !== 'string' || !/^https?:\/\//i.test(url)) continue;
          const link = document.createElement('a'); link.href = url; link.textContent = url; link.target = '_blank'; link.rel = 'noreferrer';
          const p = document.createElement('p'); p.appendChild(link); $('sources').appendChild(p);
        }
        $('evidence').textContent = JSON.stringify(bundle.enzyme_evidence, null, 2);
      } catch (err) { error(err, selectHypothesis, token); }
    }
    function filterHypotheses() {
      const mode = $('enzymeFilter').value;
      const rows = (targetData?.hypotheses || []).filter(h => mode === 'all' || (mode === 'candidate') === h.has_candidate_enzyme_evidence);
      options($('hypothesisSelect'), rows.map((h, i) => [h.id, `${i + 1}. ${h.source_reaction_ids.join(' / ')} · ${h.has_candidate_enzyme_evidence ? 'candidate evidence' : 'enzyme gap'}`]));
      $('hypothesisSelect').disabled = !rows.length;
      $('hypothesisCount').textContent = `${rows.length} of ${targetData?.hypotheses.length || 0} hypotheses`;
      selectHypothesis();
    }
    async function selectTarget() {
      const token = ++generation; clearGraph(); targetData = null; $('retry').hidden = true;
      $('hypothesisSelect').disabled = true; options($('hypothesisSelect'), []); $('hypothesisCount').textContent = '';
      const target = index.targets.find(t => t.cannabisdb_id === $('targetSelect').value);
      if (!target) { $('selectionTitle').textContent = 'No matching compound'; $('targetState').textContent = ''; $('message').textContent = 'Try a different name or CannabisDB ID.'; return; }
      $('selectionTitle').textContent = `${target.label} · ${target.cannabisdb_id}`;
      $('targetState').textContent = `${target.status}; ${target.structure_status}. ${target.next_step}`;
      if (!target.hypothesis_count) { $('message').textContent = 'No balanced net-production hypothesis for this exact target. The gap is retained; this does not establish biological absence.'; return; }
      $('message').textContent = 'Loading target hypotheses…';
      try {
        const data = await loader.target(target.cannabisdb_id);
        if (token !== generation) return;
        targetData = data; filterHypotheses();
      } catch (err) { error(err, selectTarget, token); }
    }
    function search(preferred) {
      const matches = matchingTargets(index.targets, $('targetSearch').value);
      const shown = matches.slice(0, 100);
      options($('targetSelect'), shown.map(t => [t.cannabisdb_id, `${t.label} · ${t.cannabisdb_id} · ${t.hypothesis_count} hypotheses`]));
      $('targetSelect').disabled = !shown.length;
      if (preferred && shown.some(t => t.cannabisdb_id === preferred)) $('targetSelect').value = preferred;
      $('searchCount').textContent = `${matches.length} matching records${matches.length > 100 ? '; showing first 100—refine search' : ''}`;
      selectTarget();
    }
    $('targetSearch').addEventListener('input', () => { if (index) search(); });
    $('targetSelect').addEventListener('change', selectTarget);
    $('enzymeFilter').addEventListener('change', () => { if (targetData) filterHypotheses(); });
    $('hypothesisSelect').addEventListener('change', selectHypothesis);
    async function start() {
      const token = ++generation;
      try {
        index = await loader.index();
        if (token !== generation) return;
        const s = index.summary;
        $('metrics').textContent = `${s.carbon_bearing_target_status_counts['net-production-hypotheses-found']} / ${s.carbon_bearing_target_records} carbon-bearing records have one-step hypotheses · ${s.cannabisdb_records} total records retained · not CO₂ pathway completeness`;
        if (s.carbon_bearing_targets_with_candidate_enzyme_evidence !== undefined) $('metrics').textContent += ` · ${s.carbon_bearing_targets_with_candidate_enzyme_evidence} carbon-bearing targets have candidate enzyme evidence (not confirmed activity)`;
        const requested = new URLSearchParams(location.search).get('target');
        const initial = index.targets.find(t => t.cannabisdb_id === requested) || index.targets.find(t => t.label.toLowerCase() === 'eugenol' && t.hypothesis_count) || index.targets.find(t => t.hypothesis_count);
        $('targetSearch').value = initial?.cannabisdb_id || ''; search(initial?.cannabisdb_id);
      } catch (err) { error(err, start, token); }
    }
    start();
  }
  const api = {matchingTargets, project, createLoader, mount};
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.HypothesisView = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
