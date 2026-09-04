"""Turn unresolved pathway records into explicit, testable hypotheses."""

from __future__ import annotations

import json
import re
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
            "reaction_class": reaction.get("attributes", {}).get("reactionClass"),
            "conditions": reaction.get("attributes", {}).get("conditions"),
            "reaction_smiles": reaction.get("attributes", {}).get("reactionSmiles"),
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


def _non_enzymatic_test_plan(reaction: dict) -> list[dict]:
    inputs = ", ".join(p.get("label") or p["compound_id"] for p in reaction.get("reactants", []))
    products = ", ".join(p.get("label") or p["compound_id"] for p in reaction.get("products", []))
    return [
        {"step": "controlled_decarboxylation", "method": "Incubate the isolated acidic cannabinoid under a temperature/time matrix with matched unheated controls and oxygen-controlled conditions.", "readout": f"LC-MS/MS quantification of {inputs} loss and {products} formation against authentic standards.", "claim_boundary": "This validates chemical conversion under assay conditions, not endogenous plant-pathway flux."},
        {"step": "released_carbon_validation", "method": "Quantify evolved CO2 or use a labeled carboxyl-carbon standard in a closed-vessel time course; retain the RDKit atom assignment as a structural hypothesis.", "readout": "Stoichiometric CO2 release coupled to neutral-cannabinoid formation and the expected carbon count.", "claim_boundary": "Stable-isotope evidence is required to confirm atom origin; structural mapping alone is not isotope tracing."},
        {"step": "plant_context", "method": "Compare acid/neutral ratios in fresh, dried, and deliberately heated Cannabis tissue with extraction blanks and recovery controls.", "readout": "Condition-dependent conversion consistent with the isolated-compound assay, while separating post-harvest chemistry from in-planta production."},
    ]


def _target_test_plan(target: dict) -> list[dict]:
    return [
        {"step": "identity_validation", "method": "Verify the CannabisDB structure against an authentic standard or orthogonal NMR/MS/MS evidence and resolve stereochemistry, protonation, and tautomer state.", "readout": "An exact identity match to a Terpedia/ChEBI entity, or a documented unresolved identity with alternatives retained."},
        {"step": "producer_reaction_discovery", "method": "Search Terpedia reactions, enzyme families, and the Cannabis proteome for a balanced producer reaction for the exact target structure.", "readout": "A source-linked balanced reaction with all substrates, products, and candidate proteins represented."},
        {"step": "13CO2_lineage_validation", "method": "Pulse-label Cannabis tissue with 13CO2 and measure isotopologue enrichment of the target by LC-MS/MS; use compartment and time-course controls.", "readout": "Carbon incorporation consistent with the proposed CO2 lineage; this tests carbon origin but does not by itself identify the enzyme."},
    ]


_SPECIALTY_NAME = re.compile(r"cannab|tetrahydrocannabin|cannabidiol|cannabiger|cannabichrom|cannabinol|cannabicycl|cannabielso|cannabifuran|cannabitriol|cannabid|cannabivarin|cannabistilbene|cannabisativine", re.IGNORECASE)


def _target_hypotheses(lineage_path: Path | None, compounds_path: Path | None) -> list[dict]:
    if not lineage_path or not lineage_path.exists():
        return []
    lineage = json.loads(lineage_path.read_text())
    compounds = {}
    if compounds_path and compounds_path.exists():
        compounds = {c["id"]: c for c in json.loads(compounds_path.read_text()).get("compounds", [])}
    targets = []
    for target in lineage.get("targets", []):
        reason = target.get("reason") or "unresolved-target-status"
        blockers = [] if target.get("status") == "supported" else [reason]
        compound = compounds.get(target["cannabisdb_id"], {})
        names = [compound.get("label", ""), *compound.get("aliases", [])]
        specialty = any(_SPECIALTY_NAME.search(name) for name in names)
        review_priority = "high" if specialty else "medium" if target.get("identity_status") == "exact" else "normal"
        priority_reason = "explicit Cannabis/cannabinoid specialty name" if specialty else "exact structural identity" if target.get("identity_status") == "exact" else "catalog target requiring identity resolution"
        targets.append({
            "hypothesis_id": f"target:{target['cannabisdb_id']}",
            "hypothesis_type": "metabolite-identity-and-co2-lineage",
            "status": target.get("status", "unresolved"),
            "cannabisdb_id": target["cannabisdb_id"],
            "terpedia_id": target.get("terpedia_id"),
            "identity_status": target.get("identity_status"),
            "label": compound.get("label", target["cannabisdb_id"]),
            "aliases": compound.get("aliases", []),
            "formula": compound.get("formula"),
            "smiles": compound.get("smiles"),
            "source_url": compound.get("source_url"),
            "review_priority": review_priority,
            "priority_reason": priority_reason,
            "carbon_atom_count": target.get("carbon_atom_count"),
            "reachable_carbon_atoms": target.get("reachable_carbon_atoms", 0),
            "unresolved_carbon_atoms": max((target.get("carbon_atom_count") or 0) - (target.get("reachable_carbon_atoms") or 0), 0),
            "entity_product_carbon_atoms": target.get("entity_product_carbon_atoms"),
            "reversible_upper_bound_reachable_carbon_atoms": target.get("reversible_upper_bound_reachable_carbon_atoms", 0),
            "carbon_lineage_blocker": reason,
            "claim": f"Test whether {target['cannabisdb_id']} is an endogenous Cannabis metabolite and whether its carbon atoms are traceable to inorganic CO2 through a documented pathway.",
            "blocking_causes": blockers,
            "proposed_tests": _target_test_plan(target),
            "source": str(lineage_path),
            "claim_boundary": "This is a target-validation hypothesis. Catalog presence is not evidence of endogenous biosynthesis, and CO2 isotope incorporation does not by itself establish a particular enzyme or reaction route.",
        })
    return targets


def build_test_hypotheses(queue_path: Path, network_path: Path, output: Path, lineage_path: Path | None = None, compounds_path: Path | None = None) -> dict:
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
    for reaction in reactions.values():
        if reaction.get("reaction_class") != "non-enzymatic-decarboxylation":
            continue
        hypotheses.append({
            "hypothesis_id": f"test:non-enzymatic:{reaction['id']}",
            "hypothesis_type": "non-enzymatic-conversion-validation",
            "status": "candidate",
            "reaction_id": reaction["id"],
            "reaction": reaction,
            "claim": f"The source-linked acidic cannabinoid conversion represented by {reaction['id']} occurs under defined chemical conditions.",
            "candidate_proteins": [],
            "blocking_causes": ["Cannabis-relevant condition dependence and in-planta occurrence require experimental validation"],
            "proposed_tests": _non_enzymatic_test_plan(reaction),
            "source": reaction.get("source_url"),
            "queue_kind": "non-enzymatic-reaction",
            "claim_boundary": "This is a chemical-conversion validation hypothesis. It is not an enzyme hypothesis and does not establish endogenous biosynthesis or pathway flux.",
        })
    targets = _target_hypotheses(lineage_path, compounds_path)
    report = {"schema": "cannabis-carbon.testable-hypotheses.v2", "source_queue": str(queue_path), "source_network": str(network_path), "source_lineage": str(lineage_path) if lineage_path else None, "summary": {"total": len(hypotheses), "candidate": sum(h["status"] == "candidate" for h in hypotheses), "blocked": sum(h["status"] == "blocked" for h in hypotheses), "with_reaction": sum(bool(h["reaction_id"]) for h in hypotheses), "with_candidate_proteins": sum(bool(h["candidate_proteins"]) for h in hypotheses), "non_enzymatic_conversion_hypotheses": sum(h["hypothesis_type"] == "non-enzymatic-conversion-validation" for h in hypotheses), "metabolite_targets_total": len(targets), "metabolite_targets_supported": sum(h["status"] == "supported" for h in targets), "metabolite_targets_candidate": sum(h["status"] == "candidate" for h in targets), "metabolite_targets_unresolved": sum(h["status"] == "unresolved" for h in targets)}, "hypotheses": hypotheses, "metabolite_target_hypotheses": targets, "claim_boundary": "The report turns unresolved graph records, source-linked non-enzymatic conversions, and every catalog target into falsifiable experimental hypotheses. It does not promote any candidate to a confirmed enzyme or pathway."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    return report["summary"]
