# Exact upstream gaps in four acetate routes

`phase1-acetate-precursor-gaps.json` audits the full 18,221-equation alcohol-acetate model, retaining all original directions, exact identities and external species. Its four target obstructions were independently replayed with exact fractions against every allowed directed step. Each nonnegative weighted pool includes the target and its alcohol precursor; no allowed reaction increases that weighted pool. Hence positive net target production cannot occur without internal-pool depletion in this pinned model. This is not a minimal cut, unique missing enzyme or proof of biological absence.

| Target record | Exact alcohol supply gap | Size of target obstruction |
| --- | --- | ---: |
| CDB000192 | One NADH-dependent producing step exists, but its stereospecified cyclohexanone precursor is unsupplied. The weighted family also includes alternate alcohol and ring-expanded/opened forms. | 7 compounds |
| CDB000283 | The only producing step reduces p-ethylbenzaldehyde while consuming peroxide and producing oxygen. The aldehyde/alcohol/acetate pool has no net source. This direction is a model assumption, not established plant chemistry. | 3 compounds |
| CDB000321 | The exact cannabispirol alcohol has no producing step. Adding its acetylation leaves the unsupplied alcohol/acetate pair. | 2 compounds |
| CDB000585 | NADH-dependent alcohol production exists, but the exact S-configured citronellal precursor is unsupplied. Its aldehyde, alcohol, acid and acetate form the weighted family. | 4 compounds |

No new equation or coverage gain is claimed. Detailed source reactions and all required substrates are retained in the report. The skill's all-input requirement prevents an existing alcohol-producing edge from being mistaken for a complete route.

## Literature lead for the S-citronellal gap

Martinelli and colleagues, *Plant Physiology* 194 (2024), 1006–1023, published online 13 October 2023, [DOI 10.1093/plphys/kiad550](https://doi.org/10.1093/plphys/kiad550), report three pelargonium PRISE-family citral reductases. Their abstract describes in-vitro characterization, predominantly S-citronellal or racemic products, differing NADH/NADPH preferences, and RNAi support for a role in citronellol biosynthesis. This supports a substrate-specific research lead in another plant, not Cannabis activity or enantiopure S production.

The abstract and introduction were inspected on 2026-09-06. Exact enzyme-by-substrate-isomer and cofactor assignments still need the results and methods. Subsequent targeted page retrievals and repository PDF retrieval were unsuccessful; do not claim a complete methods review. No citral reduction edge has been added from this lead. Next: resolve geranial versus neral, exact cofactor forms and measured product distributions before instantiating balanced reactions. Preserve minor stereoisomer products and organism boundaries; do not generalize an enzyme-family annotation into Cannabis substrate activity.

Follow-up: direct retrieval from the publisher's minimal-article endpoint recovered the relevant results, Figure 7 caption and enzymatic-assay methods. `citral-reductase-evidence.json` records the cofactor-specific product evidence; `phase1-citral-reductase-review.json` binds it to three distinct model identities. The unresolved point is now individual citral-isomer attribution, not availability of those methods. No single-isomer edge or enantiopure-output claim was added.
