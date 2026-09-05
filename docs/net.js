(function(root) {
  'use strict';
  function project(bundle, targetId, stepId = '') {
    const target = bundle.targets.find(t => t.cannabisdb_id === targetId);
    if (!target) throw new Error('Unknown target record');
    const certificate = bundle.certificates.find(c => c.compound_id === target.certificate_compound_id);
    if (!certificate) return {target, certificate: null, nodes: [], edges: [], steps: []};
    const compounds = new Map(bundle.compounds.map(c => [c.id, c]));
    const reactions = new Map(bundle.reactions.map(r => [r.id, r]));
    const evidence = new Map(bundle.enzyme_evidence.map(e => [e.id, e]));
    const selected = certificate.steps.filter(s => !stepId || s.step_id === stepId);
    if (!selected.length) throw new Error('Selected reaction is not in this certificate');
    const cids = new Set(), edges = [], steps = [];
    for (const step of selected) {
      const reaction = reactions.get(step.reaction_id);
      if (!reaction) throw new Error('Missing full reaction');
      if (!['hypothetical-left-to-right', 'hypothetical-right-to-left'].includes(step.direction_mode)) throw new Error('Unknown direction mode');
      const forward = step.direction_mode === 'hypothetical-left-to-right';
      const inputs = reaction[forward ? 'left' : 'right'], outputs = reaction[forward ? 'right' : 'left'];
      const attached = reaction.enzyme_evidence_ids.map(id => {if (!evidence.has(id)) throw new Error('Missing evidence'); return evidence.get(id);});
      const proteins = [...new Set(attached.flatMap(e => [...(e.screened_proteins || []).map(p => p.accession), ...(e.enzyme_ids || []), ...(e.enzyme_evidence?.screened_homology_proteins || []), ...(e.enzyme_evidence?.direction_unresolved_family_proteins || [])]))];
      steps.push({...step, reaction, required_inputs: inputs, outputs, evidence: attached});
      for (const m of [...inputs, ...outputs]) {if (!(m.coefficient > 0)) throw new Error('Invalid coefficient'); cids.add(m.compound_id);}
      inputs.forEach((a, i) => outputs.forEach((b, j) => edges.push({data: {
        id: `${step.step_id}:${i}:${j}`, source: a.compound_id, target: b.compound_id,
        kind: 'reaction-projection', step_id: step.step_id, reaction_id: step.reaction_id,
        direction_mode: step.direction_mode, extent: step.extent, required_inputs: inputs, outputs,
        enzyme_evidence_ids: reaction.enzyme_evidence_ids, candidate_protein_ids: proteins,
        is_completion_sensitivity: !!reaction.is_completion_sensitivity,
        missing_candidate_evidence: !!reaction.missing_candidate_evidence,
        claim_boundary: 'Every input is required. Projected arrows are not separate reactions, atom flow or a startup sequence.'
      }})));
    }
    const nodes = [...cids].sort().map(id => {
      const c = compounds.get(id); if (!c) throw new Error('Missing required compound');
      const pool = certificate.zero_net_internal_participants.includes(id);
      return {data: {id, label: c.labels?.[0] || `${c.formula} · ${id.split(':')[1]?.slice(0, 6) || id}`,
        compound: c, is_target: id === target.compound_id, is_pool: pool,
        role: certificate.external_net_consumption[id] ? 'input' : certificate.net_exports[id] ? 'output' : pool ? 'pool' : 'internal',
        net_consumption: certificate.external_net_consumption[id] || '0', net_export: certificate.net_exports[id] || '0'}};
    });
    return {target, certificate, nodes, edges, steps};
  }
  function matchingTargets(targets, query, scope) {
    const q = query.trim().toLowerCase();
    return targets.filter(t => (scope === 'all' || (t.certificate_compound_id && (scope !== 'enzyme-gaps' || t.missing_candidate_reaction_ids?.length))) && `${t.label} ${t.cannabisdb_id}`.toLowerCase().includes(q));
  }
  function createLoader(fetcher, folder = 'net-view') {
    if (!['net-view', 'completion-net-view', 'catalog-net-view'].includes(folder)) throw new Error('Invalid scenario folder');
    return async function() {
      const response = await fetcher(`data/${folder}/index.json`, {cache: 'no-cache'});
      if (!response.ok) throw new Error(`Manifest unavailable (HTTP ${response.status})`);
      const manifest = await response.json();
      if (manifest.file !== 'bundle.json' || !/^[a-f0-9]{64}$/.test(manifest.sha256)) throw new Error('Invalid bundle manifest');
      const data = await fetcher(`data/${folder}/bundle.json?v=${manifest.sha256.slice(0, 16)}`);
      if (!data.ok) throw new Error(`Net-conversion data unavailable (HTTP ${data.status})`);
      return data.json();
    };
  }
  function mount() {
    const scenario = new URLSearchParams(location.search).get('scenario');
    const folder = scenario === 'catalog' ? 'catalog-net-view' : scenario === 'completions' ? 'completion-net-view' : 'net-view';
    const $ = id => document.getElementById(id), loader = createLoader((...args) => fetch(...args), folder);
    if (typeof cytoscape !== 'function') { $('netMessage').textContent = 'The graph library could not load. Reload the page or use the downloadable certificates below.'; return; }
    let bundle, current, generation = 0;
    const cy = cytoscape({container: $('netCy'), elements: [], layout: {name: 'preset'}, style: [
      {selector: 'node', style: {'label': 'data(label)', 'background-color': '#aebdd2', 'color': '#e8eef8', 'font-size': 16, 'text-wrap': 'wrap', 'text-max-width': 120, 'text-valign': 'bottom', 'text-margin-y': 7, 'width': 40, 'height': 40}},
      {selector: 'node[role="input"]', style: {'background-color': '#77b9ef'}},
      {selector: 'node[role="output"]', style: {'background-color': '#65d6a0'}},
      {selector: 'node[role="pool"]', style: {'background-color': '#f2bd65'}},
      {selector: 'node[?is_target]', style: {'border-width': 5, 'border-color': '#fff'}},
      {selector: 'edge', style: {'line-color': '#c78cff', 'target-arrow-color': '#c78cff', 'target-arrow-shape': 'triangle', 'curve-style': 'bezier', 'line-style': 'dashed', 'width': 2}},
      {selector: 'edge[?is_completion_sensitivity]', style: {'line-color': '#f2a65a', 'target-arrow-color': '#f2a65a', 'width': 4}},
      {selector: 'edge[?missing_candidate_evidence]', style: {'line-color': '#ff7777', 'target-arrow-color': '#ff7777', 'width': 4}},
      {selector: '.muted', style: {'opacity': .25}}
    ]});
    function options(select, rows) {select.replaceChildren(...rows.map(([value, label]) => {const e = document.createElement('option'); e.value = value; e.textContent = label; return e;}));}
    function lines(id, title, values) {const h = document.createElement('h2'); h.textContent = title; const ul = document.createElement('ul'); for (const text of values) {const li = document.createElement('li'); li.textContent = text; ul.appendChild(li);} $(id).replaceChildren(h, ul);}
    function clear() {cy.elements().remove(); current = null; for (const id of ['netCounts','netBalance','netEquation','netSources','netEvidence','netDetails']) $(id).textContent = '';}
    function label(cid) {const c = bundle.compounds.find(c => c.id === cid); return c?.labels?.[0] || `${c?.formula || ''} (${cid})`;}
    function describe(step) {
      const equation = side => side.map(m => `${m.coefficient} × ${label(m.compound_id)}`).join(' + ');
      lines('netEquation', 'Full directed equation', [equation(step.required_inputs) + ' → ' + equation(step.outputs), `Relative extent: ${step.extent}`, 'Hypothetical direction; all listed inputs are required.']);
      $('netEvidence').textContent = JSON.stringify(step.evidence.length ? step.evidence : {
        status: 'No candidate enzyme evidence in this snapshot', reaction_id: step.reaction_id,
        claim_boundary: 'Chemistry-only step. Determine whether it is enzymatic, spontaneous or a catalog transformation; do not infer activity from connectivity.'
      }, null, 2);
      const heading = document.createElement('h2'); heading.textContent = 'Reaction sources'; $('netSources').replaceChildren(heading);
      for (const url of [...new Set(step.reaction.sources.flatMap(s => s.source_urls || []))]) {
        if (!/^https?:\/\//i.test(url)) continue;
        const p = document.createElement('p'), a = document.createElement('a'); a.href = url; a.textContent = url; a.target = '_blank'; a.rel = 'noopener noreferrer'; p.appendChild(a); $('netSources').appendChild(p);
      }
    }
    function highlight() {cy.elements().removeClass('muted'); if ($('poolHighlight').value === 'pools') cy.nodes().filter(n => !n.data('is_pool')).addClass('muted'); else if ($('poolHighlight').value === 'enzyme-gaps') cy.edges().filter(e => !e.data('missing_candidate_evidence')).addClass('muted');}
    function draw() {
      clear(); if (!bundle || !$('netTarget').value) return;
      current = project(bundle, $('netTarget').value, $('netReaction').value);
      $('netTitle').textContent = `${current.target.label} · ${current.target.cannabisdb_id}`;
      $('netStatus').textContent = `Net result: ${current.target.net_status}. Zero-pool startup: ${current.target.startup_status}.`;
      if (!current.certificate) {$('netMessage').textContent = 'No net-conversion certificate for this record. This is a model/evidence gap, not proof of biological absence.'; return;}
      const cert = current.certificate;
      $('netMessage').textContent = bundle.view_boundary || 'Exact net balance; physiological pathway and startup remain unestablished.';
      const details = ['Net inputs: ' + Object.entries(cert.external_net_consumption).map(([c,n]) => `${n} × ${label(c)}`).join(' + '),
        'Net products: ' + Object.entries(cert.net_exports).map(([c,n]) => `${n} × ${label(c)}`).join(' + '),
        'Zero-net internal participants (pool origin unresolved): ' + cert.zero_net_internal_participants.map(label).join(', ')];
      lines('netBalance', 'Complete certificate balance', details);
      cy.add([...current.nodes, ...current.edges]); cy.layout({name:'cose', animate:false, fit:true, padding:45, nodeRepulsion:20000, idealEdgeLength:100, numIter:600}).run(); highlight();
      $('netCounts').textContent = `${current.nodes.length} compounds · ${current.steps.length} directed reactions shown of ${cert.steps.length} · ${current.edges.length} projected arrows (not additional reactions)`;
      if (current.steps.length === 1) describe(current.steps[0]); else $('netEquation').textContent = 'Select an edge or a reaction above to inspect its complete equation and enzyme evidence.';
    }
    function selectTarget() {
      const target = bundle?.targets.find(t => t.cannabisdb_id === $('netTarget').value);
      const cert = bundle?.certificates.find(c => c.compound_id === target?.certificate_compound_id);
      options($('netReaction'), [['', 'Whole net conversion'], ...(cert?.steps || []).map((s,i) => [s.step_id, `${i+1}. ${bundle.reactions.find(r=>r.id===s.reaction_id).sources[0]?.source_reaction_id || s.reaction_id} · extent ${s.extent}`])]);
      $('netReaction').disabled = !cert; $('netReaction').value = ''; draw();
    }
    function search(preferred) {
      if (!bundle) return;
      const rows = matchingTargets(bundle.targets, $('netSearch').value, $('netScope').value), shown = rows.slice(0,100);
      const selected = preferred || $('netTarget').value;
      if (selected && rows.some(t=>t.cannabisdb_id===selected) && !shown.some(t=>t.cannabisdb_id===selected)) shown.push(rows.find(t=>t.cannabisdb_id===selected));
      options($('netTarget'), shown.map(t=>[t.cannabisdb_id,`${t.label} · ${t.cannabisdb_id}`]));
      $('netTarget').disabled = !shown.length; $('netTarget').value = shown.some(t=>t.cannabisdb_id===selected) ? selected : shown[0]?.cannabisdb_id || '';
      $('netMatches').textContent = `${rows.length} matches; ${shown.length} shown. Refine the search to find any record.`;
      if (!shown.length) {clear(); $('netTitle').textContent = 'No matching target'; $('netStatus').textContent = ''; $('netMessage').textContent = 'No records match this search and scope.'; options($('netReaction'), []); $('netReaction').disabled = true;} else selectTarget();
    }
    async function load() {
      const token = ++generation; clear(); bundle = null; $('netRetry').hidden = true; $('netTarget').disabled = true; $('netReaction').disabled = true;
      $('netTitle').textContent = 'Loading net-conversion evidence…'; $('netStatus').textContent = ''; $('netMessage').textContent = '';
      try {
        const loaded = await loader(); if (token !== generation) return; bundle = loaded;
        if ($('netBoundary') && bundle.view_boundary) $('netBoundary').textContent = bundle.view_boundary + ' ' + bundle.claim_boundary;
        const evidenceLabel = bundle.view_scenario === 'full-catalog-chemistry-only' ? 'chemistry-only net certificates (enzyme gaps included)' : 'candidate-linked net certificates';
        $('netMetrics').textContent = `${bundle.summary.target_status_counts['exact-net-conversion-hypothesis']} / ${bundle.summary.target_records} target records have ${evidenceLabel} · not confirmed pathway completeness`;
        const requested = new URLSearchParams(location.search).get('target');
        if (requested) {
          $('netScope').value = 'all';
          if (!bundle.targets.some(t => t.cannabisdb_id === requested)) $('netSearch').value = requested;
        }
        search(requested || bundle.targets.find(t=>t.label==='Limonene' && t.certificate_compound_id)?.cannabisdb_id);
      } catch(error) {if(token!==generation)return; clear(); $('netMessage').textContent = error.message; $('netRetry').hidden = false;}
    }
    cy.on('tap','node,edge',event=>{const data=event.target.data(); $('netDetails').textContent=JSON.stringify(data,null,2); if(data.step_id){const step=current?.steps.find(s=>s.step_id===data.step_id);if(step)describe(step);}});
    $('netFit').addEventListener('click',()=>cy.fit(cy.elements(),45)); $('poolHighlight').addEventListener('change',highlight);
    $('netSearch').addEventListener('input',()=>search()); $('netScope').addEventListener('change',()=>search());
    $('netTarget').addEventListener('change',selectTarget); $('netReaction').addEventListener('change',draw); $('netRetry').addEventListener('click',load);
    load(); return {load};
  }
  root.NetView = {project, matchingTargets, createLoader, mount};
})(typeof window !== 'undefined' ? window : globalThis);
