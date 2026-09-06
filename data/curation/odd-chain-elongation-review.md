# Odd-chain elongation sensitivity

The Terpedia C18→C20 acyl-CoA elongation sequence supplies four exact reference equations: malonyl-CoA condensation, NADPH reduction, dehydration and NADPH reduction. `phase1-odd-chain-elongation.json` retains the complete source records and proposes C17→C19 and C19→C21 homolog cycles. Only the terminal saturated chain length changes. CoA stereochemistry, specified intermediate stereochemistry, unsaturation, exact cofactors, charges and coefficients are preserved. All eight proposed equations pass independent element/isotope/charge balance checks.

These are chain-length hypotheses, not target-substrate assays, Cannabis enzyme assignments or established physiological directions. Newly added steps are forward-only assumptions. The full four-step net consumes one starter acyl-CoA, one malonyl-CoA, two NADPH and three protons; it produces one elongated acyl-CoA, CoA, CO₂, water and two NADP species in the exact charged forms retained in the report. The three intermediate acyl-CoAs cancel.

## Full-inventory result

`phase1-odd-chain-net.json` tests all 6,220 original records with unchanged identities and the same permissive 102-species external boundary, with CO₂ its sole carbon-containing input. It adds eight equations and six exact compound structures: 18,211 balanced equations and 33,532 allowed directed steps. All 2,729 prior structural certificates replay exactly. **No additional certificate or covered record is obtained.** Counts remain 2,732 conditional covered records, 2,993 missing-producer records, 460 solver-reported infeasible records and 35 exchange-species inventory records; 3,219 records participate in balanced equations. These categories are distinct from confirmed Cannabis biological coverage.

## Remaining supply obstruction

The previous model's obstruction weights no longer apply after elongation edges are added. New exact nonnegative-weight proofs in `phase1-odd-chain-obstructions.json` check every allowed directed step independently:

- Octadecane (CDB000155): an 18-compound weighted set includes C17 starter acyl-CoA and related C17 species, the C19 elongation intermediates and downstream C18 alkane.
- Eicosane (CDB000157): a 22-compound weighted set additionally includes the C21 elongation intermediates and downstream C20 alkane.

Each allowed step has nonpositive weighted net change. Thus these target pools cannot be produced net without depleting another internal pool in this pinned model. Adding downstream elongation has extended the dependency but has not supplied the upstream odd-chain pool. This is a model obstruction, not proof of biological absence or a unique missing enzyme.

Next: inspect exact C17 precursor-producing chemistry and its shorter-chain inputs, including carrier identity and stereochemistry, before proposing additional elongation or initiation reactions. Do not import odd-chain organic precursors as an unexplained shortcut. The separately running restricted-medium calculation uses the older selenium-forward model; none of these results establish a minimum medium, startup synthesis, growth, photon/energy balance or compartmental transport.

Artifacts are reproducible locally. GCP upload/read-back verification and a future Pages release are separate publication steps; generating static JSON does not establish deployment.
