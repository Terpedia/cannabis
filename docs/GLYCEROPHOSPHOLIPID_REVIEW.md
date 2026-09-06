# PG/PGP reaction completeness sensitivity

This batch adds exact-structure, balanced reaction hypotheses, not confirmed Cannabis pathways. Atom tracing and enzyme assignments are not claimed.

## Identity and reaction evidence

The full 6,220-record inventory was screened against pinned generic Rhea product scaffolds. Neutral-to-source-charge proposals match 79 PG and 78 PGP records without changing heavy atoms, bonds, isotopes or encoded stereochemistry. Of these, 151 previously lacked producing equations and six already had conditional net certificates. Each proton exchange is a separate identity-preserving reaction hypothesis with explicit H+ stoichiometry. Generic scaffold provenance is not an exact ChEBI mapping, measured pKa, tissue-pH assignment or curated protonation reaction.

The source-forward equations are RHEA:33752 (PGP hydrolysis), RHEA:12594 (PGP synthesis), RHEA:16230 (CDP-DAG synthesis), RHEA:19710 (sn-2 acylation), and RHEA:15326 (sn-1 acylation). All required water, phosphate, CTP, CMP, diphosphate, glycerol phosphate, acyl-CoA, CoA and H+ participants are retained wherever present in the source equation. No organic precursor is introduced as an external seed.

Recursive expansion examines 423 structures and produces 559 balanced proposals: 157 proton exchanges, 79 PGP hydrolyses, 102 PGP syntheses, 102 CDP-DAG syntheses, 102 sn-2 acylations and 17 sn-1 acylations. All 157 source-charge target forms have producing hypotheses. The 21 terminal precursor structures are present in the parent inventory; presence alone was not counted as supply. Independent forward bond-edit tests reconstruct every proposed synthesis/hydrolysis product, including exact stereochemistry.

## Full-network result

Merging the proposals adds 293 distinct equations and joins 266 existing equations, for 18,082 balanced equations. Existing joins retain their inherited direction constraints; new synthesis/hydrolysis reactions are forward-only and proton-exchange reversibility is an explicit sensitivity assumption.

| Metric | Parent | This sensitivity |
| --- | ---: | ---: |
| Inventory records | 6,220 | 6,220 |
| Conditional CO2 net-conversion records | 2,510 | 2,638 |
| No producing equation | 3,248 | 3,097 |
| Solver-reported infeasible | 427 | 450 |
| Explicit exchange species | 35 | 35 |
| Balanced reaction participants | 2,964 | 3,115 |

The 128 newly certified records are 51 PG and 77 PGP structures. All 2,510 parent certificates are retained. Another 23 PG records acquire producing reactions but remain infeasible under the current model; these are not coverage gains. Each new certificate passes exact balance and direction replay under the unchanged exchange boundary, with CO2 the sole external carbon input.

Certificates permit regenerated pre-existing pools. They do not establish zero-pool startup, energetics, compartment transport, physiological flux or Cannabis-specific enzyme activity. Published graph coverage should only be updated after the matching graph artifacts are tested and deployed.

## Remaining gaps and tests

- Diagnose all required inputs for the 23 newly connected but infeasible PG targets; do not assume acyl-CoA or other organic donor supply.
- Resolve CDP-DAG source identity and stereochemistry independently before proposing identity bridges. No CDP-DAG inventory records matched this exact charge-state audit.
- Test candidate PG/PGP substrate specificity experimentally using chemically resolved lipid standards. Generic Rhea scaffolds do not establish acceptance of every acyl-chain combination.
- Review protonation assumptions against compartment-specific conditions before interpreting model routes physiologically.

Machine-readable reports, exact source snapshot hashes and complete reaction records are in `data/reports/phase1-glycerophospholipid-{speciation,synthesis,net}.json`. GCP exports use immutable versioned tables in `terpedia-489015.terpedia_core`; verification receipts record full-row readback checks separately from scientific validation.
