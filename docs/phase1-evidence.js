/* Exact reaction-ID/SMARTS joins; no reaction-name or family fallback. */
(function (root) {
  'use strict';
  const key = row => JSON.stringify([row.reaction_id, row.reaction_smarts]);
  function makeIndex(rows) {
    const index = new Map();
    for (const row of rows) {
      const k = key(row);
      if (index.has(k)) throw new Error('Duplicate Phase 1 variant');
      index.set(k, row);
    }
    return index;
  }
  function select(edges, index, status) {
    const edgeIds = new Set(), nodeIds = new Set(), variants = new Set();
    for (const edge of edges) {
      const record = index.get(key(edge));
      if (!record || (status !== 'all' && record.evidence_status !== status)) continue;
      edgeIds.add(edge.id); nodeIds.add(edge.source); nodeIds.add(edge.target);
      variants.add(key(record));
    }
    return {edgeIds, nodeIds, variants};
  }
  function mount(cy) {
    const control = document.getElementById('phase1Status');
    const metrics = document.getElementById('phase1Metrics');
    const matches = document.getElementById('phase1Matches');
    let index = null, pending = false;
    function apply() {
      if (!index) return;
      const selection = select(cy.edges().map(e => e.data()), index, control.value);
      const active = control.value !== 'all';
      cy.batch(() => {
        cy.edges().forEach(e => e.toggleClass('phase1-dim', active && !selection.edgeIds.has(e.id())));
        cy.nodes().forEach(n => n.toggleClass('phase1-dim', active && !selection.nodeIds.has(n.id())));
      });
      matches.textContent = `${selection.variants.size} matching report variants on loaded graph edges. Load the three-hop candidate layer to include all expansion edges. Other filters still apply.`;
    }
    cy.style().selector('.phase1-dim').style('opacity', 0.1).update();
    control.addEventListener('change', apply);
    cy.on('add', () => {
      if (!index || pending) return;
      pending = true;
      queueMicrotask(() => { pending = false; apply(); });
    });
    cy.on('tap', 'edge', event => {
      const record = index?.get(key(event.target.data()));
      if (!record) return;
      const p = document.createElement('p');
      p.textContent = `Phase 1 variant: ${record.evidence_status}; balance: ${record.balance_status}; sequence search: ${record.search_status}. Screened proteins: ${record.screened_proteins.join(', ') || 'none'}. Core enzyme annotations: ${record.core_enzyme_ids.join(', ') || 'none'}. These are not confirmed Cannabis activities or pathways.`;
      document.getElementById('details').appendChild(p);
    });
    fetch('data/phase1-map-evidence.json').then(response => {
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      return response.json();
    }).then(report => {
      index = makeIndex(report.rows);
      const s = report.summary, counts = s.balance_status_counts;
      metrics.textContent = `Candidate expansion: ${s.balanced_variants}/${s.reaction_variants} variants balanced; ${counts.imbalanced || 0} imbalanced; ${counts.not_auditable || 0} not auditable. ${s.balanced_variants_with_candidate_enzyme_evidence}/${s.balanced_variants} balanced variants have candidate enzyme evidence; ${s.balanced_variants_without_candidate_enzyme_evidence} lack it. Not whole-metabolome or confirmed-pathway completeness.`;
      control.disabled = false;
      apply();
    }).catch(error => {
      metrics.textContent = 'Phase 1 evidence could not load. The base map remains available; reload to retry.';
      console.error(error);
    });
  }
  const api = {key, makeIndex, select, mount};
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.Phase1Evidence = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
