const {test}=require('node:test'), assert=require('node:assert/strict'), fs=require('node:fs'), vm=require('node:vm');
const script=fs.readFileSync(__dirname+'/../docs/completions.js','utf8'), context={}; vm.runInNewContext(script,context);
const {project,matching,createLoader}=context.CompletionView;
const bundle=JSON.parse(fs.readFileSync(__dirname+'/../docs/data/completion-view/bundle.json'));
test('all target completion graphs retain full coefficients and explicit inferred participants',()=>{
  let count=0;
  for(const t of bundle.targets) for(const id of t.completion_ids){
    const view=project(bundle,t.cannabisdb_id,id), h=view.hypothesis;
    const expected=new Set([...h.left,...h.right].map(p=>p.compound_id));
    assert.deepEqual(new Set(view.nodes.map(n=>n.data.id)),expected);
    assert.equal(view.edges.length,view.inputs.length*view.outputs.length);
    for(const e of view.edges){assert.equal(e.data.required_inputs,view.inputs);assert.equal(e.data.outputs,view.outputs);assert.equal(e.data.hypothesis_id,id);}
    for(const n of view.nodes) assert.equal(n.data.inferred,h.inferred_inorganic_participants_in_MARTS_orientation.flat().some(p=>p.compound_id===n.data.id));
    count++;
  }
  assert.ok(count>=67); assert.equal(matching(bundle,'','all').length,6220);
  assert.equal(matching(bundle,'','hypotheses').length,67);
  assert.equal(matching(bundle,'','hypotheses','with').length,62);
  const without=matching(bundle,'','hypotheses','without');
  assert.equal(without.length,8); // Three targets have both supported and unsupported alternatives.
  assert.equal(without.filter(t=>t.completion_ids.every(id=>!bundle.completions.find(h=>h.id===id).protein_evidence.has_candidate_lead)).length,5);
  const gap=bundle.targets.find(t=>!t.completion_ids.length);assert.equal(project(bundle,gap.cannabisdb_id).nodes.length,0);
  assert.throws(()=>project(bundle,gap.cannabisdb_id,bundle.completions[0].id),/belong/);
  const archived=bundle.completions.find(h=>h.protein_evidence.archived_source_screen?.has_archived_candidate_lead);
  assert.ok(archived); assert.equal(archived.protein_evidence.has_candidate_lead,true);
  assert.ok(archived.protein_evidence.representative_alignments.some(a=>a.reference_accession.startsWith('UPI')));
  const archivedTarget=bundle.targets.find(t=>t.completion_ids.includes(archived.id));
  assert.ok(matching(bundle,archivedTarget.cannabisdb_id,'hypotheses','with').length);
});
test('loader revalidates and versions local data and rejects external paths',async()=>{
  const calls=[];let fail=true;
  const loader=createLoader(async(url,options)=>{calls.push({url,options});if(fail){fail=false;return {ok:false};}return {ok:true,json:async()=>url.endsWith('index.json')?{file:'bundle.json',sha256:'a'.repeat(64)}:bundle};});
  await assert.rejects(loader(),/unavailable/);assert.equal(await loader(),bundle);
  assert.equal(calls[1].options.cache,'no-cache');assert.ok(calls[2].url.endsWith('?v='+'a'.repeat(16)));
  await assert.rejects(createLoader(async()=>({ok:true,json:async()=>({file:'https://invalid/file',sha256:'a'.repeat(64)})}))(),/manifest/);
});
test('controls select a target, clear gaps, and retry without unsafe HTML',async()=>{
  class Field{constructor(){this.value='';this.children=[];this.textContent='';}addEventListener(e,f){this[e]=f;}replaceChildren(){this.children=[];this.value='';}append(x){this.children.push(x);if(!this.value&&x.value)this.value=x.value;}set innerHTML(x){throw Error('unsafe HTML');}}
  const names=['Search','Scope','Target','Choice','Metrics','Matches','Title','Message','Equation','Status','Counts','Evidence','Details','Sources','Cy','Fit','Retry','Protein','ProteinFilter'];
  const fields=Object.fromEntries(names.map(n=>['completion'+n,new Field()]));fields.completionScope.value='hypotheses';
  fields.completionProteinFilter.value='all';
  let failure=true, destroyed=0;
  const env={URLSearchParams,location:{search:''},document:{getElementById:id=>fields[id],createElement:()=>new Field()},
    cytoscape:()=>({on(){},fit(){},destroy(){destroyed++;}}),fetch:async url=>{if(failure){failure=false;return {ok:false};}return {ok:true,json:async()=>url.endsWith('index.json')?{file:'bundle.json',sha256:'a'.repeat(64)}:bundle};}};
  vm.runInNewContext(script,env);await env.CompletionView.mount();assert.equal(fields.completionRetry.hidden,false);
  await fields.completionRetry.click();assert.ok(fields.completionEquation.textContent.includes('[inferred]'));
  fields.completionProteinFilter.value='without';fields.completionProteinFilter.change();
  assert.equal(JSON.parse(fields.completionEvidence.textContent).protein_evidence.has_candidate_lead,false);
  fields.completionProteinFilter.value='with';fields.completionProteinFilter.change();
  assert.equal(JSON.parse(fields.completionEvidence.textContent).protein_evidence.has_candidate_lead,true);
  fields.completionProteinFilter.value='all';fields.completionProteinFilter.change();
  const archived=bundle.completions.find(h=>h.protein_evidence.archived_source_screen?.has_archived_candidate_lead);
  fields.completionSearch.value=bundle.targets.find(t=>t.completion_ids.includes(archived.id)).cannabisdb_id;fields.completionSearch.input();
  fields.completionChoice.value=archived.id;fields.completionChoice.change();
  assert.ok(fields.completionProtein.children.some(c=>c.textContent.includes('Archived sequence identity is not functional annotation')));
  assert.ok(JSON.parse(fields.completionEvidence.textContent).protein_evidence.archived_source_screen.has_archived_candidate_lead);
  const gap=bundle.targets.find(t=>!t.completion_ids.length);fields.completionScope.value='all';fields.completionSearch.value=gap.cannabisdb_id;fields.completionSearch.input();
  assert.equal(fields.completionEquation.textContent,'');assert.ok(fields.completionStatus.textContent.startsWith('No completion'));assert.ok(destroyed>0);
});
