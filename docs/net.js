(function(root) {
  'use strict';
  function applyEvidence(base, supplement) {
    if(base.view_scenario !== 'full-catalog-chemistry-only' || supplement.schema !== 'cannabis-carbon.phase1-catalog-evidence.v1') throw Error('Invalid evidence scenario');
    const existing=new Set(base.enzyme_evidence.map(e=>e.id)), added=new Map();
    for(const e of supplement.enzyme_evidence) {
      const r=base.reactions.find(r=>r.id===e.reaction_id);
      if(!r || r.enzyme_evidence_ids.length || existing.has(e.id) || added.has(e.reaction_id) || !e.screened_proteins.length) throw Error('Invalid or duplicate evidence reaction');
      existing.add(e.id); added.set(e.reaction_id,e);
    }
    const certUpdates=new Map(supplement.certificate_updates.map(u=>[u.compound_id,u]));
    const targetUpdates=new Map(supplement.target_updates.map(u=>[u.cannabisdb_id,u]));
    if(certUpdates.size!==supplement.certificate_updates.length || targetUpdates.size!==supplement.target_updates.length) throw Error('Duplicate evidence update');
    const update = (item, change) => {
      const before=item.missing_candidate_reaction_ids, after=before.filter(rid=>!added.has(rid));
      if(after.length===before.length) {if(change) throw Error('Unchanged record in supplement');return item;}
      if(!change || change.compound_id!==item.compound_id || JSON.stringify(change.baseline_missing_candidate_reaction_ids)!==JSON.stringify(before) || JSON.stringify(change.missing_candidate_reaction_ids)!==JSON.stringify(after)) throw Error('Evidence gap update mismatch');
      return {...item,baseline_missing_candidate_reaction_ids:before,missing_candidate_reaction_ids:after,
        evidence_class:after.length?'chemistry-only-with-enzyme-gaps':'candidate-linked-net-hypothesis',has_new_catalog_evidence:true};
    };
    const certificates=base.certificates.map(c=>update(c,certUpdates.get(c.compound_id)));
    const targets=base.targets.map(t=>update(t,targetUpdates.get(t.cannabisdb_id)));
    if(certificates.filter(c=>c.has_new_catalog_evidence).length!==certUpdates.size || targets.filter(t=>t.has_new_catalog_evidence).length!==targetUpdates.size) throw Error('Unknown evidence update');
    return {...base,certificates,targets,
      reactions:base.reactions.map(r=>added.has(r.id)?{...r,enzyme_evidence_ids:[added.get(r.id).id],missing_candidate_evidence:false,is_new_catalog_candidate:true}:r),
      enzyme_evidence:[...base.enzyme_evidence,...supplement.enzyme_evidence],
      evidence_summary:supplement.summary,view_boundary:supplement.view_boundary,
      claim_boundary:base.claim_boundary+' '+supplement.claim_boundary};
  }
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
      const directionAudit = reaction.certificate_direction_annotations?.[step.direction_mode] || null;
      const attached = reaction.enzyme_evidence_ids.map(id => {if (!evidence.has(id)) throw new Error('Missing evidence'); return evidence.get(id);});
      const proteins = [...new Set(attached.flatMap(e => [...(e.screened_proteins || []).map(p => p.accession), ...(e.enzyme_ids || []), ...(e.enzyme_evidence?.screened_homology_proteins || []), ...(e.enzyme_evidence?.direction_unresolved_family_proteins || [])]))];
      steps.push({...step, reaction, required_inputs: inputs, outputs, evidence: attached});
      for (const m of [...inputs, ...outputs]) {if (!(m.coefficient > 0)) throw new Error('Invalid coefficient'); cids.add(m.compound_id);}
      inputs.forEach((a, i) => outputs.forEach((b, j) => edges.push({data: {
        id: `${step.step_id}:${i}:${j}`, source: a.compound_id, target: b.compound_id,
        kind: 'reaction-projection', step_id: step.step_id, reaction_id: step.reaction_id,
        direction_mode: step.direction_mode, extent: step.extent, required_inputs: inputs, outputs,
        hypothesis_type: reaction.hypothesis_type || null,
        source_reaction_id: reaction.source_reaction_id || null,
        source_url: reaction.source_url || null,
        source_mapping_evidence: reaction.source_mapping_evidence || [],
        speciation_type: reaction.speciation_type || null,
        direction_status: reaction.direction_status || null,
        balance_status: reaction.balance_status || null,
        light_requirement_annotation: reaction.light_requirement_annotation || null,
        certificate_direction_annotation: directionAudit,
        enzyme_evidence_ids: reaction.enzyme_evidence_ids, candidate_protein_ids: proteins,
        is_completion_sensitivity: !!reaction.is_completion_sensitivity,
        is_route_sensitivity: !!reaction.is_route_sensitivity,
        hypothesis_assumptions: reaction.hypothesis_assumptions || [],
        missing_candidate_evidence: !!reaction.missing_candidate_evidence,
        is_new_catalog_candidate: !!reaction.is_new_catalog_candidate,
        has_unreviewed_reference: attached.some(e => e.evidence_class?.includes('unreviewed') || e.reference_matches?.some(r=>r.review_status==='unreviewed')),
        direction_review_id: reaction.direction_review?.id || (directionAudit?.review_flags.length ? step.step_id : null),
        direction_review_warning: reaction.direction_review?.warning || null,
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
    return targets.filter(t => (scope === 'all' || (t.certificate_compound_id && (scope !== 'enzyme-gaps' || (t.missing_candidate_reaction_count ?? t.missing_candidate_reaction_ids?.length)))) && `${t.label} ${t.cannabisdb_id}`.toLowerCase().includes(q));
  }
  function shardedView(base, fetcher, folder) {
    if (base.schema !== 'cannabis-carbon.sharded-net-view.v1' || !Array.isArray(base.targets) || base.targets.length !== base.summary?.target_records || !base.certificate_files) throw Error('Invalid sharded index');
    const validRef = (r, pattern) => r && pattern.test(r.file) && /^[a-f0-9]{64}$/.test(r.sha256) && Number.isSafeInteger(r.bytes) && r.bytes > 0;
    if (!validRef(base.shared_chemistry, /^chemistry\.json$/)) throw Error('Invalid shared chemistry reference');
    const targets = new Map(base.targets.map(t=>[t.cannabisdb_id,t]));
    if(targets.size !== base.targets.length) throw Error('Duplicate target accession');
    for (const [cid,ref] of Object.entries(base.certificate_files)) if (!/^structure:[a-f0-9]{64}$/.test(cid) || !validRef(ref,/^certificates\/[a-f0-9]{64}\.json$/)) throw Error('Invalid certificate reference');
    for (const t of base.targets) {
      if (!Number.isSafeInteger(t.missing_candidate_reaction_count) || t.missing_candidate_reaction_count < 0) throw Error('Invalid enzyme-gap count');
      if (t.certificate_compound_id && (t.certificate_compound_id !== t.compound_id || !base.certificate_files[t.compound_id])) throw Error('Invalid target certificate identity');
    }
    async function read(ref) {
      const response = await fetcher(`data/${folder}/${ref.file}?v=${ref.sha256.slice(0,16)}`);
      if (!response.ok) throw Error(`Pathway data unavailable (HTTP ${response.status})`);
      const bytes = await response.arrayBuffer();
      if (bytes.byteLength !== ref.bytes) throw Error('Pathway data size mismatch');
      const digest = Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes)),n=>n.toString(16).padStart(2,'0')).join('');
      if (digest !== ref.sha256) throw Error('Pathway data checksum mismatch');
      return JSON.parse(new TextDecoder().decode(bytes));
    }
    let shared;
    const certificates = new Map();
    const state = {...base, certificates:[], reactions:[], compounds:[], enzyme_evidence:[], loadTarget};
    async function loadTarget(id) {
      const target = targets.get(id);
      if (!target) throw Error('Unknown target record');
      const cid = target.certificate_compound_id;
      if (!cid) return {...state};
      if (!shared) shared = read(base.shared_chemistry).then(data=>{
        for (const key of ['reactions','compounds','enzyme_evidence']) if (!Array.isArray(data[key])) throw Error('Invalid shared chemistry');
        return data;
      }).catch(error=>{shared=null;throw error;});
      if (!certificates.has(cid)) {
        const promise=read(base.certificate_files[cid]).then(cert=>{
          if(cert.compound_id!==cid || !Array.isArray(cert.steps)) throw Error('Certificate identity mismatch');
          return cert;
        }).catch(error=>{certificates.delete(cid);throw error;});
        certificates.set(cid,promise);
        if(certificates.size>16) certificates.delete(certificates.keys().next().value);
      }
      const [chemistry,cert] = await Promise.all([shared,certificates.get(cid)]);
      const byReaction = new Map(chemistry.reactions.map(r=>[r.id,r]));
      const missing = new Set();
      for (const step of cert.steps) {
        const reaction = byReaction.get(step.reaction_id);
        if (!reaction || base.forbidden_step_ids?.includes(step.step_id)) throw Error('Missing or forbidden certificate reaction');
        if (!reaction.enzyme_evidence_ids.length) missing.add(step.reaction_id);
      }
      if (missing.size!==target.missing_candidate_reaction_count) throw Error('Certificate enzyme-gap count mismatch');
      return {...state,...chemistry,certificates:[cert],targets:base.targets.map(t=>t.compound_id===cid ? {...t,missing_candidate_reaction_ids:[...missing].sort()} : t)};
    }
    return state;
  }
  function createLoader(fetcher, folder = 'net-view', comparison = 'CHI-and-FNSII') {
    if(['local-speciation-net-view','selenium-net-view','ketone-stereo-net-view','c17-net-view'].includes(folder)) return async()=>{
      const response=await fetcher(`data/${folder}/index.json`,{cache:'no-cache'});
      if(!response.ok) throw Error(`Manifest unavailable (HTTP ${response.status})`);
      const manifest=await response.json();
      if(manifest.file!=='bundle.json'|| !/^[a-f0-9]{64}$/.test(manifest.sha256)) throw Error('Invalid bundle manifest');
      const data=await fetcher(`data/${folder}/bundle.json?v=${manifest.sha256.slice(0,16)}`);
      if(!data.ok) throw Error(`Net-conversion data unavailable (HTTP ${data.status})`);
      return shardedView(await data.json(),fetcher,folder);
    };
if (!['pg-named-net-view', 'glycerophospholipid-net-view', 'amino-phospholipid-net-view', 'source-mapped-protonation-net-view', 'cardiolipin-net-view', 'glycerolipid-precursors-net-view', 'phosphatidate-hydrolysis-net-view', 'triglyceride-symmetry-net-view', 'triglyceride-net-view', 'chemistry-net-view', 'fnsii-net-view', 'net-view', 'completion-net-view', 'catalog-net-view', 'expanded-net-view', 'purine-net-view', 'purine-restricted-net-view', 'thiolase-net-view', 'thiolase-restricted-net-view', 'remaining-net-view', 'remaining-restricted-net-view'].includes(folder)) throw new Error('Invalid scenario folder');
    if (folder === 'fnsii-net-view' && !['baseline', 'CHI-only', 'FNSII-only', 'CHI-and-FNSII'].includes(comparison)) throw new Error('Invalid sensitivity comparison');
    const sourceFolder = folder === 'remaining-restricted-net-view' ? 'remaining-net-view' : folder === 'thiolase-restricted-net-view' ? 'thiolase-net-view' : folder === 'purine-restricted-net-view' ? 'purine-net-view' : folder;
    return async function() {
      const response = await fetcher(`data/${sourceFolder}/index.json`, {cache: 'no-cache'});
      if (!response.ok) throw new Error(`Manifest unavailable (HTTP ${response.status})`);
      const manifest = await response.json();
      if (manifest.file !== 'bundle.json' || !/^[a-f0-9]{64}$/.test(manifest.sha256)) throw new Error('Invalid bundle manifest');
      const data = await fetcher(`data/${sourceFolder}/bundle.json?v=${manifest.sha256.slice(0, 16)}`);
      if (!data.ok) throw new Error(`Net-conversion data unavailable (HTTP ${data.status})`);
      const base = await data.json();
      if (folder === 'pg-named-net-view') {
        if (!['original_result','alternative_baseline_result','alternative_extended_result'].includes(comparison)) throw Error('Invalid identity comparison');
        const selected = base.scenario_options?.find(s => s.id === comparison);
        if (base.view_scenario !== 'pg-paired-identity' || !selected || selected.targets?.length !== 23 || !Array.isArray(selected.certificates)) throw Error('Invalid paired identity scenario');
        return {...base,...selected,view_boundary:selected.title + ': ' + base.view_boundary};
      }
      if (folder === 'fnsii-net-view') {
        const selected = base.scenario_options?.find(s => s.id === comparison);
        if (base.view_scenario !== 'fnsii-route-sensitivity' || base.model_eligible !== false || !selected || !Array.isArray(selected.certificates) || !Array.isArray(selected.targets)) throw Error('Invalid sensitivity scenario');
        return {...base, ...selected, view_boundary: `${comparison}: ${base.view_boundary}`};
      }
      if (folder === 'remaining-restricted-net-view') {
        if(base.view_scenario !== 'remaining-candidates' || base.restricted_scenario?.id !== 'eight-reverse-steps-forbidden') throw Error('Invalid restricted scenario');
        return {...base, ...base.restricted_scenario, view_boundary:base.restricted_boundary};
      }
      if (folder === 'thiolase-restricted-net-view') {
        if(base.view_scenario !== 'thiolase-candidates' || base.restricted_scenario?.id !== 'eight-reverse-steps-forbidden') throw Error('Invalid restricted scenario');
        return {...base, ...base.restricted_scenario, view_boundary:base.restricted_boundary};
      }
      if (folder === 'purine-restricted-net-view') {
        if(base.view_scenario !== 'purine-candidates' || base.restricted_scenario?.id !== 'five-reverse-steps-forbidden') throw Error('Invalid restricted scenario');
        return {...base, ...base.restricted_scenario, view_boundary:base.restricted_boundary};
      }
      if (!manifest.evidence) return base;
      const extra = manifest.evidence;
      if(folder !== 'catalog-net-view' || extra.file !== 'evidence.json' || !/^[a-f0-9]{64}$/.test(extra.sha256)) throw Error('Invalid evidence manifest');
      const responseExtra = await fetcher(`data/${folder}/evidence.json?v=${extra.sha256.slice(0,16)}`);
      if(!responseExtra.ok) throw Error(`Evidence supplement unavailable (HTTP ${responseExtra.status})`);
      const supplement = await responseExtra.json();
      if(supplement.source_sha256?.['data/reports/phase1-catalog-net-gaps.json'] !== manifest.sha256) throw Error('Evidence snapshot mismatch');
      return applyEvidence(base,supplement);
    };
  }
  function mount() {
    const scenario = new URLSearchParams(location.search).get('scenario');
    const folder = scenario === 'hydrolysis' ? 'phosphatidate-hydrolysis-net-view' : scenario === 'symmetry' ? 'triglyceride-symmetry-net-view' : scenario === 'triglycerides' ? 'triglyceride-net-view' : scenario === 'chemistry' ? 'chemistry-net-view' : scenario === 'fnsii' ? 'fnsii-net-view' : scenario === 'remaining' ? 'remaining-net-view' : scenario === 'remaining-restricted' ? 'remaining-restricted-net-view' : scenario === 'thiolase' ? 'thiolase-net-view' : scenario === 'thiolase-restricted' ? 'thiolase-restricted-net-view' : scenario === 'purine' ? 'purine-net-view' : scenario === 'purine-restricted' ? 'purine-restricted-net-view' : scenario === 'expanded' ? 'expanded-net-view' : scenario === 'catalog' ? 'catalog-net-view' : scenario === 'completions' ? 'completion-net-view' : 'net-view';
const selectedFolder = scenario === 'ketone' ? 'ketone-stereo-net-view' : scenario === 'selenium' ? 'selenium-net-view' : scenario === 'speciation' ? 'local-speciation-net-view' : scenario === 'pg-named' ? 'pg-named-net-view' : scenario === 'glycerophospholipids' ? 'glycerophospholipid-net-view' : scenario === 'aminos' ? 'amino-phospholipid-net-view' : scenario === 'protonation' ? 'source-mapped-protonation-net-view' : scenario === 'cardiolipins' ? 'cardiolipin-net-view' : scenario === 'precursors' ? 'glycerolipid-precursors-net-view' : folder;
    const $ = id => document.getElementById(id), loader = createLoader((...args) => fetch(...args), scenario === 'c17' ? 'c17-net-view' : selectedFolder, new URLSearchParams(location.search).get('comparison') || (scenario === 'pg-named' ? 'alternative_extended_result' : 'CHI-and-FNSII'));
    if (typeof cytoscape !== 'function') { $('netMessage').textContent = 'The graph library could not load. Reload the page or use the downloadable certificates below.'; return; }
    let bundle, current, generation = 0, selectionGeneration = 0;
    const cy = cytoscape({container: $('netCy'), elements: [], layout: {name: 'preset'}, style: [
      {selector: 'node', style: {'label': 'data(label)', 'background-color': '#aebdd2', 'color': '#e8eef8', 'font-size': 16, 'text-wrap': 'wrap', 'text-max-width': 120, 'text-valign': 'bottom', 'text-margin-y': 7, 'width': 40, 'height': 40}},
      {selector: 'node[role="input"]', style: {'background-color': '#77b9ef'}},
      {selector: 'node[role="output"]', style: {'background-color': '#65d6a0'}},
      {selector: 'node[role="pool"]', style: {'background-color': '#f2bd65'}},
      {selector: 'node[?is_target]', style: {'border-width': 5, 'border-color': '#fff'}},
      {selector: 'edge', style: {'line-color': '#c78cff', 'target-arrow-color': '#c78cff', 'target-arrow-shape': 'triangle', 'curve-style': 'bezier', 'line-style': 'dashed', 'width': 2}},
      {selector: 'edge[?is_completion_sensitivity]', style: {'line-color': '#f2a65a', 'target-arrow-color': '#f2a65a', 'width': 4}},
      {selector: 'edge[?missing_candidate_evidence]', style: {'line-color': '#ff7777', 'target-arrow-color': '#ff7777', 'width': 4}},
      {selector: 'edge[?is_new_catalog_candidate]', style: {'line-color': '#65c8ff', 'target-arrow-color': '#65c8ff', 'width': 4}},
      {selector: '.muted', style: {'opacity': .25}}
    ]});
    function options(select, rows) {select.replaceChildren(...rows.map(([value, label]) => {const e = document.createElement('option'); e.value = value; e.textContent = label; return e;}));}
    function lines(id, title, values) {const h = document.createElement('h2'); h.textContent = title; const ul = document.createElement('ul'); for (const text of values) {const li = document.createElement('li'); li.textContent = text; ul.appendChild(li);} $(id).replaceChildren(h, ul);}
    function clear() {cy.elements().remove(); current = null; for (const id of ['netCounts','netBalance','netEquation','netSources','netEvidence','netDetails']) $(id).textContent = '';}
    function label(cid) {const c = bundle.compounds.find(c => c.id === cid); return c?.labels?.[0] || `${c?.formula || ''} (${cid})`;}
    function describe(step) {
      const equation = side => side.map(m => `${m.coefficient} × ${label(m.compound_id)}`).join(' + ');
      lines('netEquation', 'Full directed equation', [equation(step.required_inputs) + ' → ' + equation(step.outputs), `Relative extent: ${step.extent}`, 'Hypothetical direction; all listed inputs are required.', ...(step.reaction.hypothesis_assumptions || []), ...(step.reaction.direction_review ? [step.reaction.direction_review.warning, ...step.reaction.direction_review.discriminating_tests] : [])]);
      const audit = step.reaction.certificate_direction_annotations?.[step.direction_mode];
      if (audit) {
        const warnings = {
          'peroxide-consuming-oxygen-producing-direction-review': 'This direction consumes peroxide and produces oxygen. Review its redox and physiological feasibility; balance alone does not establish it.',
          'co2-consuming-step-not-automatically-photosynthetic-fixation': 'This step consumes CO₂; that alone does not establish photosynthetic fixation.'
        };
        const p = document.createElement('p');
        p.textContent = [audit.evidence_class === 'proposed-gap-filling-reaction' ? 'Proposed gap-filling chemistry.' : 'Catalog-linked chemistry; not proof of Cannabis activity.', ...audit.review_flags.map(f=>warnings[f] || f), 'This audit does not establish a Cannabis enzyme or physiological direction.'].join(' ');
        $('netEquation').appendChild(p);
      }
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
    function highlight() {cy.elements().removeClass('muted'); if ($('poolHighlight').value === 'pools') cy.nodes().filter(n => !n.data('is_pool')).addClass('muted'); else if ($('poolHighlight').value === 'route-sensitivity') cy.edges().filter(e => !e.data('is_route_sensitivity')).addClass('muted'); else if ($('poolHighlight').value === 'enzyme-gaps') cy.edges().filter(e => !e.data('missing_candidate_evidence')).addClass('muted'); else if ($('poolHighlight').value === 'direction-review') cy.edges().filter(e => !e.data('direction_review_id')).addClass('muted'); else if ($('poolHighlight').value === 'unreviewed-references') cy.edges().filter(e => !e.data('has_unreviewed_reference')).addClass('muted');}
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
      const selection = ++selectionGeneration, loadGeneration = generation;
      if(bundle?.loadTarget) {
        const id = $('netTarget').value;
        clear(); $('netReaction').disabled = true; $('netRetry').hidden = true;
        $('netTitle').textContent = 'Loading selected pathway…'; $('netStatus').textContent = '';
        $('netMessage').textContent = 'Loading full equations and the selected certificate…';
        return bundle.loadTarget(id).then(loaded=>{
          if(selection!==selectionGeneration || loadGeneration!==generation) return;
          bundle=loaded; renderTarget();
        }).catch(error=>{
          if(selection!==selectionGeneration || loadGeneration!==generation) return;
          clear(); $('netMessage').textContent=error.message; $('netRetry').hidden=false;
        });
      }
      renderTarget();
    }
    function renderTarget() {
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
      if (!shown.length) {++selectionGeneration; clear(); $('netTitle').textContent = 'No matching target'; $('netStatus').textContent = ''; $('netMessage').textContent = 'No records match this search and scope.'; options($('netReaction'), []); $('netReaction').disabled = true;} else return selectTarget();
    }
    async function load() {
      const token = ++generation; clear(); bundle = null; $('netRetry').hidden = true; $('netTarget').disabled = true; $('netReaction').disabled = true;
      $('netTitle').textContent = 'Loading net-conversion evidence…'; $('netStatus').textContent = ''; $('netMessage').textContent = '';
      try {
        const loaded = await loader(); if (token !== generation) return; bundle = loaded;
        if ($('netBoundary') && bundle.view_boundary) $('netBoundary').textContent = bundle.view_boundary + ' ' + bundle.claim_boundary;
        const evidenceLabel = bundle.view_scenario === 'pg-paired-identity' ? 'conditional identity-comparison certificates (no historical coverage gain)' : bundle.view_scenario === 'fnsii-route-sensitivity' ? 'conditional sensitivity certificates (assumed enzyme steps)' : ['full-catalog-chemistry-only', 'reaction-first-chemistry'].includes(bundle.view_scenario) ? 'chemistry-only net certificates (enzyme gaps included)' : 'candidate-linked net certificates';
        $('netMetrics').textContent = `${bundle.summary.target_status_counts['exact-net-conversion-hypothesis'] ?? 0} / ${bundle.summary.target_records} target records have ${evidenceLabel} · not confirmed pathway completeness`;
        if(bundle.evidence_summary) $('netMetrics').textContent += ` · ${bundle.evidence_summary.selected_certificate_targets_with_candidates_for_all_steps} selected target certificates have candidates for all steps · ${bundle.evidence_summary.remaining_missing_candidate_equations} reaction gaps remain`;
        const requested = new URLSearchParams(location.search).get('target');
        if (['fnsii-net-view', 'pg-named-net-view'].includes(selectedFolder) && !bundle.certificates.length) $('netScope').value = 'all';
        if (requested) {
          $('netScope').value = 'all';
          if (!bundle.targets.some(t => t.cannabisdb_id === requested)) $('netSearch').value = requested;
        }
        const newlyFeasible = ['expanded-candidates','purine-candidates','thiolase-candidates','remaining-candidates'].includes(bundle.view_scenario) ? bundle.targets.find(t=>t.new_net_certificate)?.cannabisdb_id : null;
        await search(requested || newlyFeasible || bundle.targets.find(t=>t.label==='Limonene' && t.certificate_compound_id)?.cannabisdb_id);
      } catch(error) {if(token!==generation)return; clear(); $('netMessage').textContent = error.message; $('netRetry').hidden = false;}
    }
    cy.on('tap','node,edge',event=>{const data=event.target.data(); $('netDetails').textContent=JSON.stringify(data,null,2); if(data.step_id){const step=current?.steps.find(s=>s.step_id===data.step_id);if(step)describe(step);}});
    $('netFit').addEventListener('click',()=>cy.fit(cy.elements(),45)); $('poolHighlight').addEventListener('change',highlight);
    $('netSearch').addEventListener('input',()=>search()); $('netScope').addEventListener('change',()=>search());
    $('netTarget').addEventListener('change',selectTarget); $('netReaction').addEventListener('change',draw); $('netRetry').addEventListener('click',load);
    load(); return {load};
  }
  root.NetView = {project, matchingTargets, createLoader, applyEvidence, mount};
})(typeof window !== 'undefined' ? window : globalThis);
