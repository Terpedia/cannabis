# Alkane hypothesis sensitivity

The full 6,220-record inventory was tested after adding eight forward-only hypotheses: four C19–C22 acyl-CoA reductions to aldehydes and four aldehyde deformylations to C18–C21 alkanes. Exact existing substrate structures and the complete NADPH/CoA/oxygen/proton/formate/water bookkeeping are retained. These are chain-length analogies, not measured Cannabis reactions or plant CER1 assignments.

Two new certificates pass exact rational replay:

- CDB000156, nonadecane: 68 reaction steps; 19 CO₂ net input per target molecule.
- CDB000158, heneicosane: 71 reaction steps; 21 CO₂ net input per target molecule.

Both use the proposed acyl-CoA reduction and deformylation. The formate coproduct is accounted for by the full network, not discarded as unexplained carbon. Both also import peroxide and iron–sulfur species under the inherited permissive boundary. They must not be represented as validated minimum-medium pathways.

CDB000155 (octadecane) and CDB000157 (eicosane) now have producing equations but remain solver-reported infeasible. The C19 acyl-CoA for octadecane exists in the catalog, showing why precursor presence alone is insufficient. The C21 acyl-CoA for eicosane is a newly proposed exact structure without an established supply route.

Summary: 18,203 balanced equations; 3,219 participating records; 2,732 conditionally covered records / 2,729 exact structural certificates; 2,993 missing-producer records; 460 solver-reported infeasible records; 35 external-species inventory records. All 2,727 prior structural certificates were replayed and retained. No exact identities, exchange permissions or inherited direction bounds were changed.

Next work: resolve odd-chain acyl-CoA supply and source-specific substrate scope; separate plant alkane chemistry from cyanobacterial deformylation; test the paired routes with curated input and energy constraints. Current IUBMB EC 1.2.1.50 describes an acyl-protein complex and separately notes acyl-CoA acceptance. It must not be used to merge protein/ACP carriers with CoA or to assign a Cannabis protein merely from an EC label.

Sources are retained in `phase1-alkane-hypotheses.json` and `phase1-alkane-precursors.json`; full target results, inherited constraints and exact certificates are in `phase1-alkane-net.json`. These results extend the ketone scenario and are separate from the in-progress medium restriction on the earlier selenium-forward model.
