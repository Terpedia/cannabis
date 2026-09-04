"""Turn unresolved pathway records into explicit, testable hypotheses."""

from __future__ import annotations

import json
from pathlib import Path


def _reaction_index(network: dict) -> dict[str, dict]:
    entities = {entity["id"]: entity for entity in network.get("entities", [])}
    statements = network.get("statements", [])
    index = {}
    for reaction in (e for e in network.get("entities", []) if e.get("type") == "biochemical_reaction"):
        rid = reaction["id"]
        rows = [s for s in statements if s.get("subjectId") == rid]
        index[rid] = {
            "id": rid,
            "label": reaction.get("label", rid),
            "equation": reaction.get("attributes", {}).get("equation"),
            "ec_numbers": reaction.get("attributes", {}).get("ecNumbers", []),
            "source_url": reaction.get("url"),
            "reactants": [{"compound_id": s.get("objectEntityId"), "coefficient": (s.get("qualifiers") or {}).get("stoichiometricCoefficient", 1), "label": entities.get(s.get("objectEntityId"), {}).get("label")} for s in rows if s.get("predicate") == "has_reactant"],
            "products": [{"compound_id": s.get("objectEntityId"), "coefficient": (s.get("qualifiers") or {}).get("stoichiometricCoefficient", 1), "label": entities.get(s.get("objectEntityId"), {}).get("label")} for s in rows if s.get("predicate") == "has_product"],
        }
    return index


def _best_hit(protein: dict) -> dict | None:
    search = protein.get("diamond_search") or protein.get("specialized_search") or {}
    hits = search.get("hits") or []
    return sorted(hits, key=lambda hit: (hit.get("evalue", float("inf")), -hit.get("bitscore", 0)))[0] if hits else None


def _test_plan(reaction_id: str | None, reaction: dict | None) -> list[dict]:
    inputs = ", ".join(f"{p['coefficient']}× {p.get('label') or p['compound_id']}" for p in (reaction or {}).get("reactants", []))
    products = ", ".join(f"{p['coefficient']}× {p.get('label') or p['compound_id']}" for p in (reaction or {}).get("products", []))
    plan = [{"step": "recombinant_assay", "method": "Express the candidate protein, purify it, and compare an active-enzyme reaction to no-enzyme and heat-inactivated controls.", "readout": f"LC-MS/MS detection of the expected reaction products: {products or 'reaction products'}.", "inputs": inputs or "Resolve the missing reaction substrates before assay design."}]
    if reaction_id == "cannabis:reaction:tks-hexanoyl-coa-to-tetraketide-coa":
        plan.append({"step": "coupled_TKS_OAC_assay", "method": "Assay candidate TKS alone and with characterized OAC using hexanoyl-CoA and malonyl-CoA.", "readout": "Detect tetraketide-CoA and/or olivetolate by LC-MS; include authentic standards where available."})
    if reaction_id == "cannabis:reaction:tetraketide-coa-to-olivetolate":
        plan.append({"step": "OAC_substrate_specificity", "method": "Provide purified tetraketide-CoA (or a validated coupled TKS preparation) to candidate OAC proteins.", "readout": "Olivetolate formation with carboxylate retention, compared with substrate-free and candidate-free controls."})
    plan.append({"step": "plant_validation", "method": "Measure candidate transcript/protein abundance in cannabinoid-producing glandular trichomes and compare chemotypes.", "readout": "Co-localization and genotype/chemotype association; this supports biological relevance but does not replace the enzyme assay."})
    return plan


def build_test_hypotheses(queue_path: Path, network_path: Path, output: Path) -> dict:
    queue = json.loads(queue_path.read_text())
    if network_path.suffix == ".gz":
        from .terpedia import load_network
        network = load_network(network_path)
    else:
        network = json.loads(network_path.read_text())
    reactions = _reaction_index(network)
    hypotheses = []
    for item in queue.get("items", []):
        reaction_id = item.get("reaction_id")
        reaction = reactions.get(reaction_id) if reaction_id else None
        candidates = []
        for protein in item.get("candidate_proteins", []):
            hit = _best_hit(protein)
            candidates.append({"protein_id": protein.get("proteinId"), "accession": protein.get("accession"), "label": protein.get("label"), "gene_symbols": protein.get("geneSymbols", []), "sequence_present": (protein.get("sequence_search") or {}).get("sequence_present"), "best_hit": hit, "candidate_origin": protein.get("candidateOrigin")})
        blockers = []
        if not candidates:
            blockers.append("no candidate protein is attached to this hypothesis")
        if reaction is None and reaction_id:
            blockers.append("reaction is not present in the working reaction index")
        if reaction and not reaction["reactants"]:
            blockers.append("reaction has no resolved reactant participants")
        hypothesis_type = "candidate-enzyme-function" if reaction_id else "reaction-or-producer-discovery"
        if item.get("kind") == "blocked-known-reaction":
            hypothesis_type = "blocked-reaction-mechanism"
        claim = f"A Cannabis protein candidate catalyzes {reaction['label']}" if reaction else "A missing or unresolved Cannabis reaction/protein hypothesis can explain the queued pathway record."
        hypotheses.append({"hypothesis_id": f"test:{item.get('id')}", "hypothesis_type": hypothesis_type, "status": "candidate" if candidates else "blocked", "reaction_id": reaction_id, "reaction": reaction, "claim": claim, "candidate_proteins": candidates, "blocking_causes": blockers, "proposed_tests": _test_plan(reaction_id, reaction), "source": item.get("source"), "queue_kind": item.get("kind"), "claim_boundary": "This is a test design and candidate ranking record. It is not evidence of enzyme activity, pathway flux, or in-vivo cannabis biosynthesis."})
    report = {"schema": "cannabis-carbon.testable-hypotheses.v1", "source_queue": str(queue_path), "source_network": str(network_path), "summary": {"total": len(hypotheses), "candidate": sum(h["status"] == "candidate" for h in hypotheses), "blocked": sum(h["status"] == "blocked" for h in hypotheses), "with_reaction": sum(bool(h["reaction_id"]) for h in hypotheses), "with_candidate_proteins": sum(bool(h["candidate_proteins"]) for h in hypotheses)}, "hypotheses": hypotheses, "claim_boundary": "The report turns unresolved graph records into falsifiable experimental hypotheses. It does not promote any candidate to a confirmed enzyme or pathway."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    return report["summary"]
