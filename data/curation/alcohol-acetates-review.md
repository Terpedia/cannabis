# Alcohol acetate reaction-completeness sensitivity

## Selection and exact identity

The full C17 scenario has 14 missing-producer target records containing an acetate ester. The structure-based pass selects neutral, unlabeled-isotope primary/secondary alcohol monoacetates only. It explicitly leaves tertiary-alcohol acetates, phenolic/enolic oxygen substrates and multiple-ester targets for separate review. These exclusions are scope boundaries, not biological impossibility claims.

Seven exact alcohol identities can be obtained by removing only the acetyl group. One, the alcohol from CDB000623 (source label `lsobornyl acetate`), is absent from the current exact-identity model and receives no reaction. Six have existing exact alcohol and acetyl-CoA inputs and receive proposed acetylation equations. Existing precursor presence does not establish net supply.

Every accepted target round-trips from its derived alcohol to the original isomeric SMILES. Protonation, specified stereochemistry and bond orders are preserved. In particular, the alcohol derived from the record labeled `BORNYL acetATE` matches the current structure labeled `(-)-Isoborneol`. This is a source-label/structure review lead, not permission to relabel or merge stereoisomers; coverage applies to the recorded structure. The distinct CDB000623 alcohol is not silently substituted.

## Reference and evidence boundary

The exact Terpedia reference is balanced equation `b2ec3259c31c217eaad784236834d1d76101c7fa619d5694b20cd6546ceb772c`, retaining RHEA:36148/36149 source records. In its canonical forward orientation, benzyl alcohol and acetyl-CoA yield benzyl acetate and CoA. Each proposal replaces only the alcohol/ester pair and preserves both exact charged CoA species and unit stoichiometry. No free acetate, ATP, water or proton is silently added to the net equation.

[IUBMB EC 2.3.1.84](https://iubmb.qmul.ac.uk/enzyme/EC2/3/1/84.html), accessed 2026-09-06, describes the general acetyl-CoA/alcohol reaction and reports a range of short-chain aliphatic substrates. This does not establish the proposed cyclic, terpenoid or aromatic alcohol scope, nor assign that EC to Cannabis proteins or to the exact benzyl reference. Full substrate and organism evidence remain separate. All six proposed equations are explicitly forward-only class analogies, not enzyme-confirmed reactions.

## Full-inventory outcome

`phase1-alcohol-acetates-net.json` retains all 2,733 prior exact structural certificates and adds two: CDB000110 (`BORNYL acetATE`, exact stored identity) and CDB000695 (isobutyl acetate). Four additional targets acquire producing reactions but remain solver-reported infeasible: CDB000192, CDB000283, CDB000321 and CDB000585. Their alcohol presence is not sufficient for a full net pathway.

Across all 6,220 records, counts are **2,738 conditionally covered records / 2,735 structural certificates**, 2,987 records without producing equations, 460 solver-reported infeasible records and 35 exchange-species inventory records. There are 18,221 balanced equations and 3,226 balanced-participating records. The increase in infeasible records reflects four missing-producer gaps moving into a distinct upstream-supply-gap category, not loss of prior coverage.

CO₂ remains the sole external carbon input; all 102 permissive exchange species and inherited direction restrictions are unchanged. Exact certificate replay does not establish initial pool synthesis, physiological uptake, light, energetics, growth or Cannabis activity. These results are not the ongoing restricted-medium calculation or a medium recipe.

The pathway-inference skill's exact-identity and all-input requirements determined the precursor checks and kept the four incomplete routes out of coverage. Next steps are to inspect those four upstream alcohol gaps, reconcile the bornyl/isobornyl source identity evidence without overwriting the original records, and test substrate-specific acetyltransfer activity. Publication to GCP and live Pages is established only by separate verification receipts.
