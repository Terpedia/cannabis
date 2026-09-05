(function(root) {
  'use strict';
  function project(bundle, targetId, completionId) {
    const target = bundle.targets.find(t => t.cannabisdb_id === targetId);
    if (!target) throw new Error('Unknown target');
    const id = completionId || target.completion_ids[0];
    if (!id) return {target, nodes: [], edges: []};
    if (!target.completion_ids.includes(id)) throw new Error('Hypothesis does not belong to target');
    const hypothesis = bundle.completions.find(h => h.id === id);
    if (!hypothesis) throw new Error('Missing completion');
    if (!['hypothetical-left-to-right', 'hypothetical-right-to-left'].includes(hypothesis.marts_forward_direction)) throw new Error('Unknown direction');
    const [inputs, outputs] = hypothesis.marts_forward_direction === 'hypothetical-left-to-right' ? [hypothesis.left, hypothesis.right] : [hypothesis.right, hypothesis.left];
    const compounds = new Map(bundle.compounds.map(c => [c.id,c]));
    const labels = new Map();
    for (const t of bundle.targets) if (!labels.has(t.compound_id)) labels.set(t.compound_id,t.label);
    labels.set(target.compound_id,target.label);
    const inferred = new Set(hypothesis.inferred_inorganic_participants_in_MARTS_orientation.flat().map(p => p.compound_id));
    const nodes = [...new Set([...inputs,...outputs].map(p => p.compound_id))].map(id => {
      const c = compounds.get(id); if (!c) throw new Error('Missing compound identity');
      return {data: {id, label: labels.get(id) || c.formula, compound: c, inferred: inferred.has(id), is_target: id === target.compound_id}};
    });
    const edges = [];
    inputs.forEach((a,i) => outputs.forEach((b,j) => edges.push({data: {
      id: `${hypothesis.id}:${i}:${j}`, source: a.compound_id, target: b.compound_id,
      hypothesis_id: hypothesis.id, reaction_id: hypothesis.balanced_equation_id,
      required_inputs: inputs, outputs, inferred_inorganic_participants: hypothesis.inferred_inorganic_participants_in_MARTS_orientation,
      direction_status: hypothesis.direction_status, evidence_class: hypothesis.status, claim_boundary: hypothesis.claim_boundary
    }})));
    const label = p => `${p.coefficient === 1 ? '' : p.coefficient + ' × '}${labels.get(p.compound_id) || compounds.get(p.compound_id).formula}${inferred.has(p.compound_id) ? ' [inferred]' : ''}`;
    return {target, hypothesis, nodes, edges, inputs, outputs, equation: inputs.map(label).join(' + ') + ' → ' + outputs.map(label).join(' + ')};
  }
  function matching(bundle, query, scope) {
    query = query.trim().toLowerCase();
    return bundle.targets.filter(t => (scope === 'all' || t.completion_ids.length) && `${t.cannabisdb_id} ${t.label}`.toLowerCase().includes(query));
  }
  function createLoader(fetcher) {
    return async function() {
      const res = await fetcher('data/completion-view/index.json', {cache:'no-cache'});
      if (!res.ok) throw new Error('Manifest unavailable');
      const manifest = await res.json();
      if (manifest.file !== 'bundle.json' || !/^[a-f0-9]{64}$/.test(manifest.sha256)) throw new Error('Invalid manifest');
      const data = await fetcher(`data/completion-view/bundle.json?v=${manifest.sha256.slice(0,16)}`);
      if (!data.ok) throw new Error('Hypothesis data unavailable');
      return data.json();
    };
  }
  function mount() {
    const doc=root.document, el=id=>doc.getElementById('completion'+id), loader=createLoader(root.fetch.bind(root));
    let bundle, cy, generation=0;
    function option(select,value,label) { const o=doc.createElement('option'); o.value=value; o.textContent=label; select.append(o); }
    function clear() {
      if(cy) {cy.destroy(); cy=null;}
      for(const name of ['Equation','Status','Counts','Evidence','Details']) el(name).textContent='';
      el('Sources').replaceChildren(); el('Title').textContent='Choose a target';
    }
    function link(parent,url,label) {
      if(!/^https?:\/\//i.test(url||'')) return;
      const a=doc.createElement('a'); a.href=url; a.textContent=label; a.target='_blank'; a.rel='noopener noreferrer';
      const p=doc.createElement('p'); p.append(a); parent.append(p);
    }
    function render() {
      clear(); if(!bundle || !el('Target').value) return;
      const view=project(bundle,el('Target').value,el('Choice').value);
      el('Title').textContent=`${view.target.label} · ${view.target.cannabisdb_id}`;
      if(!view.hypothesis) {el('Status').textContent='No completion hypothesis in this audit scope. This is not evidence of biological absence.'; return;}
      const h=view.hypothesis;
      el('Equation').textContent=view.equation;
      el('Status').textContent='Balanced stoichiometric hypothesis; enzyme, direction, product assignment and all-input supply require review.';
      el('Evidence').textContent=JSON.stringify(h,null,2);
      el('Counts').textContent=`${view.nodes.length} compounds · 1 completion hypothesis · ${view.edges.length} projected arrows`;
      const variant=bundle.variants.find(v=>v.id===h.variant_id);
      const details=doc.createElement('details'), summary=doc.createElement('summary'); summary.textContent='Original MARTS source records'; details.append(summary);
      for(const sid of variant.source_record_ids) {
        const s=bundle.MARTS_sources.find(s=>s.id===sid).source_record;
        link(details,s.source_url,s.rule_id+' · original equation remains unbalanced');
        if(s.source_uniprot_id) link(details,'https://www.uniprot.org/'+(s.source_uniprot_id.startsWith('UPI')?'uniparc/':'uniprotkb/')+encodeURIComponent(s.source_uniprot_id),s.source_uniprot_id+' · source reference only');
      }
      el('Sources').append(details);
      const refs=doc.createElement('details'), title=doc.createElement('summary'); title.textContent='Cofactor templates · no enzyme transfer'; refs.append(title);
      const urls=new Set();
      for(const t of h.reference_templates) {
        const r=bundle.reference_reactions.find(r=>r.id===t.reference_reaction_id);
        for(const s of r.sources) for(const url of s.source_urls||[]) if(!urls.has(url)) {urls.add(url); link(refs,url,s.source_reaction_id);}
      }
      el('Sources').append(refs);
      if(!root.cytoscape) {el('Message').textContent='Graph library unavailable. The full equation and source data remain available.';return;}
      cy=root.cytoscape({container:el('Cy'),elements:[...view.nodes,...view.edges],style:[
        {selector:'node',style:{label:'data(label)','font-size':16,color:'#e8eef8','text-valign':'bottom','text-margin-y':8,'background-color':'#65d6a0',width:46,height:46}},
        {selector:'node[?inferred]',style:{'background-color':'#f2bd65'}},
        {selector:'node[?is_target]',style:{'border-width':3,'border-color':'#fff'}},
        {selector:'edge',style:{'curve-style':'bezier','target-arrow-shape':'triangle','target-arrow-color':'#f2bd65','line-color':'#f2bd65','line-style':'dashed',width:2}}
      ],layout:{name:'breadthfirst',directed:true,padding:60,spacingFactor:1.4,animate:false}});
      cy.on('tap','node,edge',event=>{el('Details').textContent=JSON.stringify(event.target.data(),null,2);});
    }
    function choose() {
      el('Choice').replaceChildren();
      const target=bundle.targets.find(t=>t.cannabisdb_id===el('Target').value);
      for(const [i,id] of (target?.completion_ids||[]).entries()) option(el('Choice'),id,'Hypothesis '+(i+1)+' · inferred stoichiometry');
      el('Choice').disabled=!target?.completion_ids.length; render();
    }
    function filter(preferred) {
      const rows=matching(bundle,el('Search').value,el('Scope').value), shown=rows.slice(0,100);
      const selected=rows.find(t=>t.cannabisdb_id===preferred);
      if(selected&&!shown.includes(selected)) shown.push(selected);
      el('Target').replaceChildren(); shown.forEach(t=>option(el('Target'),t.cannabisdb_id,`${t.label} · ${t.cannabisdb_id}`));
      if(selected) el('Target').value=selected.cannabisdb_id;
      el('Target').disabled=!shown.length;
      el('Matches').textContent=`${rows.length} matches · showing ${shown.length}${rows.length>shown.length?' · refine search':''}`;
      choose();
    }
    async function load() {
      const ticket=++generation; bundle=null; clear(); el('Target').replaceChildren(); el('Choice').replaceChildren();
      el('Metrics').textContent='Loading…'; el('Matches').textContent='';
      el('Target').disabled=true; el('Choice').disabled=true; el('Retry').hidden=true; el('Message').textContent='Loading completion hypotheses…';
      try {
        const result=await loader(); if(ticket!==generation)return; bundle=result;
        el('Metrics').textContent=`${bundle.summary.completion_hypotheses} full-catalog hypotheses · ${bundle.summary.targets_with_completions} targets with hypotheses · ${bundle.targets.length} total targets`;
        el('Message').textContent='';
        const requested=new URLSearchParams(root.location?.search||'').get('target');
        if(requested) {el('Scope').value='all';el('Search').value=requested;}
        filter(requested||'CDB000078');
      } catch(error) {if(ticket!==generation)return;bundle=null;clear();el('Target').replaceChildren();el('Choice').replaceChildren();el('Target').disabled=true;el('Choice').disabled=true;el('Metrics').textContent='Data unavailable';el('Message').textContent='Unable to load hypotheses. Retry or use the report links below.';el('Retry').hidden=false;}
    }
    el('Search').addEventListener('input',()=>{if(bundle)filter();}); el('Scope').addEventListener('change',()=>{if(bundle)filter();});
    el('Target').addEventListener('change',()=>{if(bundle)choose();}); el('Choice').addEventListener('change',render);
    el('Fit').addEventListener('click',()=>{if(cy)cy.fit(undefined,50);}); el('Retry').addEventListener('click',load);
    return load();
  }
  root.CompletionView={project,matching,createLoader,mount};
})(typeof window!=='undefined'?window:globalThis);
