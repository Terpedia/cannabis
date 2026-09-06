# Carrier scaffold review

The source-carrier reconstruction contains 1,342 equation variants from 2,739
source joins. A conservation screen finds 1,305 variants with equal proposed
redox-scaffold counts, 35 with unresolved scaffold definitions, and two with
different source carrier contexts on opposite sides.

This is a conditional symbolic check, not a complete macromolecular formula
certificate. Source names supply conservation hypotheses only; carrier IDs are
never merged. The pathway-inference checks require full substrates, direction,
regeneration, and explicit startup assumptions before a pathway claim.

## Source-context exceptions

| Source reactions | One side | Other side | Required review |
| --- | --- | --- | --- |
| RHEA:80604 / 80605 | 2 Fe(II)-[cytochrome], Compound_11778 | 2 Fe(III)-[cytochrome c], Compound_14399 | Determine whether the generic carrier is intended to mean cytochrome c; preserve both IDs until source evidence resolves scope. |
| RHEA:55233 / 55234 | 2 Fe(III)-[cytochrome b561], Compound_10441 | 2 Fe(II)-[cytochrome b5], Compound_10438 | Check the source participant assignment and primary reaction evidence; these are distinct named carrier contexts. |

The table gives model-side comparisons, not physiological directions. The
source records and exact side correspondence remain in the reconstruction.
These discrepancies are not repaired by substituting free iron, equating the
carrier names, or assuming a source error. They remain unresolved.

## Consequences for the medium and pathway model

Carrier-bound metals and redox groups are internal molecular states, not free
nutrient imports. Equal symbolic scaffold counts do not establish carrier
regeneration across a complete route or show that its initial carrier pool can
be synthesized. No external carrier uptake is authorized by this screen.

Next, reconstruct the full model with the distinct source carrier species and
source-specific direction records. Preserve unresolved scaffold, polymer, and
derived-template cases separately; do not silently count their projected
small-molecule balance as full identity validation. Recompute qualified
CO2-to-target coverage only with the resulting assumptions made explicit.

Artifacts: `data/reports/phase1-carrier-scaffolds.json` and
`data/reports/phase1-carrier-reconstruction.json`. The screen covers every
restored variant but does not cover carrier contexts absent from that source
join or prove full reaction completeness for any Cannabis metabolite.
