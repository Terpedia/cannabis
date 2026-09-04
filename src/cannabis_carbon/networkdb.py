"""Build a unified, source-preserving NetworkDB snapshot."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

from rdkit import Chem

from .terpedia import load_network


def build_map_snapshot(networkdb_path: Path, output: Path, lineage_path: Path | None = None, focus_output: Path | None = None) -> dict:
    """Write a compact visualization snapshot containing reaction-connected records."""
    networkdb = json.loads(networkdb_path.read_text())
    reaction_connected_ids = {p["compound_id"] for r in networkdb.get("reactions", []) for p in r.get("reactants", []) + r.get("products", [])}
    hypothetical_connections = networkdb.get("hypothetical_connections", [])
    hypothesis_connected_ids = {h[side] for h in hypothetical_connections for side in ("substrate_compound_id", "product_compound_id")}
    lineage = json.loads(lineage_path.read_text()) if lineage_path and lineage_path.exists() else {}
    reachable_entities = set(lineage.get("reachable_carbon_entity_ids", []))
    compounds = [{**c, "reaction_connected": c.get("id") in reaction_connected_ids, "hypothesis_connected": c.get("id") in hypothesis_connected_ids, "co2_reachable": c.get("id") in reachable_entities or bool(c.get("co2_reachable_carbon_atoms"))} for c in networkdb.get("compounds", [])]
    reactions = []
    for reaction in networkdb.get("reactions", []):
        compact = {key: reaction.get(key) for key in ("id", "label", "equation", "ec_numbers", "reactants", "products", "enzyme_ids", "status", "carbon_mapping", "source_url", "directional_rhea_ids", "direction")}
        compact["candidate_proteins"] = [{"proteinId": p.get("proteinId"), "accession": p.get("accession"), "label": p.get("label")} for p in reaction.get("candidate_proteins", [])]
        reactions.append(compact)
    report = {"schema": "cannabis-carbon.network-map.v1", "source_networkdb": str(networkdb_path), "source_lineage": str(lineage_path) if lineage_path else None, "claim_boundary": "This is a compact visualization projection containing the complete compound inventory, balanced reaction records, and a separate hypothesis-edge layer. It is not proof of in-vivo flux.", "compounds": compounds, "reactions": reactions, "hypothetical_connections": hypothetical_connections, "coverage": networkdb.get("coverage", {}), "focus": {"co2_reachable_compounds": sum(c["co2_reachable"] for c in compounds), "reaction_connected_compounds": sum(c["reaction_connected"] for c in compounds), "hypothesis_connected_compounds": sum(c["hypothesis_connected"] for c in compounds), "all_inventory_compounds": len(compounds)}}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    if focus_output:
        focused_ids = {c["id"] for c in compounds if c["reaction_connected"] or c["co2_reachable"] or c["hypothesis_connected"]}
        focused_reactions = [
            reaction for reaction in reactions
            if any(p["compound_id"] in focused_ids for p in reaction.get("reactants", []) + reaction.get("products", []))
        ]
        focused = {**report, "schema": "cannabis-carbon.network-map-focus.v1", "compounds": [c for c in compounds if c["id"] in focused_ids], "reactions": focused_reactions, "hypothetical_connections": [h for h in hypothetical_connections if h["substrate_compound_id"] in focused_ids or h["product_compound_id"] in focused_ids], "focus": {**report["focus"], "focused_compounds": len(focused_ids)}}
        focus_output.parent.mkdir(parents=True, exist_ok=True)
        focus_output.write_text(json.dumps(focused, separators=(",", ":")) + "\n")
    return {"compounds": len(compounds), "reactions": len(reactions), "bytes": output.stat().st_size}


def build_networkdb(network_path: Path, compounds_path: Path, crosswalk_path: Path, output: Path, hypotheses_path: Path | None = None, genome_search_path: Path | None = None, genome_fasta_path: Path | None = None, mapping_path: Path | None = None, lineage_path: Path | None = None, atom_audit_path: Path | None = None, pubchem_path: Path | None = None, identity_set_path: Path | None = None, hypothetical_connections_path: Path | None = None) -> dict:
    network = load_network(network_path)
    directions_path = network_path.parent / "directional-reaction-overrides.json"
    directions = json.loads(directions_path.read_text()) if directions_path.exists() else {}
    catalog = json.loads(compounds_path.read_text())["compounds"]
    crosswalk = json.loads(crosswalk_path.read_text())
    pubchem = json.loads(pubchem_path.read_text()) if pubchem_path and pubchem_path.exists() else None
    pubchem_by_id = {r["cannabisdb_id"]: r for r in (pubchem or {}).get("records", [])}
    identity_set = json.loads(identity_set_path.read_text()) if identity_set_path and identity_set_path.exists() else None
    identity_set_by_id = {r["cannabisdb_id"]: r for r in (identity_set or {}).get("records", [])}
    hypothetical = json.loads(hypothetical_connections_path.read_text()) if hypothetical_connections_path and hypothetical_connections_path.exists() else None
    hypotheses = json.loads(hypotheses_path.read_text()) if hypotheses_path and hypotheses_path.exists() else {"items": []}
    hypothesis_items = hypotheses.get("items") or hypotheses.get("hypotheses", [])
    genome_search = json.loads(genome_search_path.read_text()) if genome_search_path and genome_search_path.exists() else None
    genome_by_protein = {p["proteinId"]: p for p in (genome_search or {}).get("candidate_proteins", [])}
    proteome_sequences = {}
    mapping_by_reaction = {}
    if mapping_path and mapping_path.exists():
        mapping_by_reaction = {row["reaction_id"]: row for row in json.loads(mapping_path.read_text()).get("reactions", [])}
    lineage_edges_by_reaction = {}
    lineage_edge_totals = {status: 0 for status in ("inferred", "candidate")}
    lineage_targets = {}
    if lineage_path and lineage_path.exists():
        lineage_report = json.loads(lineage_path.read_text())
        lineage_targets = {target["cannabisdb_id"]: target for target in lineage_report.get("targets", [])}
        for edge in lineage_report.get("carbon_edges", []):
            lineage_edges_by_reaction.setdefault(edge["reaction_id"], {status: 0 for status in ("inferred", "candidate")})[edge["status"]] += 1
            lineage_edge_totals[edge["status"]] += 1
    if genome_fasta_path and genome_fasta_path.exists():
        from .genome import _fasta
        proteome_sequences = _fasta(genome_fasta_path)

    def enrich_protein(protein):
        accession = protein.get("accession")
        sequence = proteome_sequences.get(accession, "")
        if sequence:
            return {**protein, "sequence_search": {"method": "reference-proteome-membership", "proteome": "UP000583929", "fasta_source": str(genome_fasta_path), "sequence_present": True, "length": len(sequence), "sha256": hashlib.sha256(sequence.encode()).hexdigest(), "claim": "Sequence presence and annotation support a candidate protein; they do not establish reaction specificity or in-vivo activity."}}
        return protein
    entities = {e["id"]: e for e in network["entities"]}
    identity_by_cdb = {row["cannabisdb"]["cannabisdb_id"]: row["terpedia_id"] for row in crosswalk["matches"]}
    candidate_identity_by_cdb = {}
    for row in crosswalk.get("candidate_matches", []):
        candidate_identity_by_cdb.setdefault(row["cannabisdb"]["cannabisdb_id"], []).append({"terpedia_id": row["terpedia_id"], "terpedia_label": row.get("terpedia_label"), "method": row.get("method"), "identity_status": "candidate"})
    proteins_by_ec = {}
    for protein in (e for e in network["entities"] if e.get("type") == "protein"):
        attrs = protein.get("attributes", {})
        record = {"proteinId": protein["id"], "accession": protein.get("identifiers", {}).get("uniprotAccession"), "label": protein.get("label"), "geneSymbols": attrs.get("geneSymbols", []), "exactEcNumbers": attrs.get("exactEcNumbers", []), "sourceUrl": protein.get("url") or f"https://www.uniprot.org/uniprotkb/{protein.get('identifiers', {}).get('uniprotAccession', '')}/entry", "candidateOrigin": "Terpedia exact EC annotation"}
        for ec in attrs.get("exactEcNumbers", []):
            proteins_by_ec.setdefault(ec, []).append(record)
    compounds = []
    for compound in catalog:
        candidates = candidate_identity_by_cdb.get(compound["id"], [])
        target = lineage_targets.get(compound["id"], {})
        pubchem_record = pubchem_by_id.get(compound["id"], {})
        compounds.append({**compound, "namespace": "cannabisdb", "identity_link": identity_by_cdb.get(compound["id"]), "identity_link_candidates": candidates, "identity_status": "exact" if compound["id"] in identity_by_cdb else "candidate" if candidates else "unresolved", "terpedia_identity_set": identity_set_by_id.get(compound["id"]), "pubchem": pubchem_record.get("pubchem"), "pubchem_status": pubchem_record.get("status", "not-queried"), "pubchem_reason": pubchem_record.get("reason"), "carbon_lineage_status": target.get("status", "unresolved"), "co2_reachable_carbon_atoms": target.get("reachable_carbon_atoms", 0), "carbon_lineage_reason": target.get("reason", "no-lineage-record")})
    for entity in network["entities"]:
        if entity.get("type") != "metabolite":
            continue
        attrs = entity.get("attributes", {})
        mol = Chem.MolFromSmiles(attrs["canonicalSmiles"]) if attrs.get("canonicalSmiles") else None
        compounds.append({"id": entity["id"], "namespace": "terpedia", "label": entity.get("label", entity["id"]), "formula": attrs.get("molecularFormula"), "smiles": attrs.get("canonicalSmiles"), "source_url": entity.get("url"), "carbon_atom_count": sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()) if mol is not None else None})
    hypothetical_connections = []
    if hypothetical:
        existing_compound_ids = {compound["id"] for compound in compounds}
        for row in hypothetical.get("connections", []):
            substrate_id = f"terpedia:identity-set:{row['normalized_substrate_terpene_id']}"
            product_id = f"terpedia:identity-set:{row['normalized_product_terpene_id']}"
            for compound_id, terpene_id, smiles, formula, inchikey, carbon_count, identity_set_key in (
                (substrate_id, row["normalized_substrate_terpene_id"], row["substrate_identity_smiles"], row["substrate_identity_formula"], row["substrate_inchikey"], row["substrate_carbon_count"], row["substrate_identity_set_key"]),
                (product_id, row["normalized_product_terpene_id"], row["product_identity_smiles"], row["product_identity_formula"], row["product_inchikey"], row["product_carbon_count"], row["product_identity_set_key"]),
            ):
                if compound_id not in existing_compound_ids:
                    compounds.append({"id": compound_id, "namespace": "terpedia_identity_set", "label": terpene_id, "formula": formula, "smiles": smiles, "inchikey": inchikey, "carbon_atom_count": carbon_count, "identity_set_key": identity_set_key, "identity_status": "identity-set-record", "source_url": "https://console.cloud.google.com/bigquery?project=terpedia-489015"})
                    existing_compound_ids.add(compound_id)
            hypothetical_connections.append({**row, "substrate_compound_id": substrate_id, "product_compound_id": product_id, "status": "candidate", "layer": "Terpedia hypothesis edge", "claim_boundary": "This source-directed connection is a testable hypothesis and is not included in balanced reaction or CO2 lineage counts."})
    candidate_by_reaction = {}
    for item in hypothesis_items:
        if item.get("reaction_id"):
            proteins = []
            for candidate in item.get("candidate_proteins", []):
                p = {**candidate, "proteinId": candidate.get("proteinId") or candidate.get("protein_id"), "geneSymbols": candidate.get("geneSymbols") or candidate.get("gene_symbols", []), "candidateOrigin": candidate.get("candidateOrigin") or candidate.get("candidate_origin")}
                protein_id = p.get("proteinId")
                proteins.append({**enrich_protein(p), "sequence_search": genome_by_protein.get(protein_id, {}).get("sequence_search") or enrich_protein(p).get("sequence_search"), "diamond_search": genome_by_protein.get(protein_id, {}).get("diamond_search")})
            candidate_by_reaction.setdefault(item["reaction_id"], []).extend(proteins)
    reactions = []
    for entity in network["entities"]:
        if entity.get("type") != "biochemical_reaction":
            continue
        reaction_id = entity["id"]
        statements = [s for s in network["statements"] if s.get("subjectId") == reaction_id]
        participants = lambda predicate: [{"compound_id": s["objectEntityId"], "coefficient": (s.get("qualifiers") or {}).get("stoichiometricCoefficient", 1), "compartment": (s.get("qualifiers") or {}).get("compartment")} for s in statements if s.get("predicate") == predicate]
        enzyme_statements = [s for s in network["statements"] if s.get("predicate") in ("catalyzes", "maps_to_reaction", "has_catalytic_activity") and s.get("objectEntityId") == reaction_id]
        enzymes = sorted({s.get("subjectId") for s in enzyme_statements})
        attrs = entity.get("attributes", {})
        candidate_proteins = candidate_by_reaction.get(reaction_id, [])
        annotation_candidates = [{**enrich_protein(p), "candidateOrigin": "Terpedia reaction EC to protein exact-EC join"} for ec in attrs.get("ecNumbers", []) for p in proteins_by_ec.get(ec, [])]
        candidate_proteins = list({p.get("proteinId"): p for p in candidate_proteins + annotation_candidates}.values())
        mapping = mapping_by_reaction.get(reaction_id)
        mapping_counts = {status: sum(m.get("status") == status for m in (mapping or {}).get("mappings", [])) for status in ("inferred", "candidate", "ambiguous", "unresolved")}
        carbon_status = "unavailable" if mapping is None else "unresolved" if mapping_counts["unresolved"] else "ambiguous" if mapping_counts["ambiguous"] else "candidate" if mapping_counts["candidate"] else "inferred"
        carbon_mapping = {"status": carbon_status, "product_carbon_atoms": (mapping or {}).get("product_carbon_atom_count"), "counts": mapping_counts, "lineage_edge_counts": lineage_edges_by_reaction.get(reaction_id, {status: 0 for status in ("inferred", "candidate")}), "source": str(mapping_path) if mapping_path else None, "lineage_source": str(lineage_path) if lineage_path else None}
        raw_reactants, raw_products = participants("has_reactant"), participants("has_product")
        direction = directions.get(reaction_id, {})
        reactants, products = (raw_products, raw_reactants) if direction.get("orientation") == "reverse_master" else (raw_reactants, raw_products)
        reactions.append({"id": reaction_id, "label": entity.get("label", reaction_id), "equation": attrs.get("equation"), "reaction_smiles": attrs.get("reactionSmiles"), "ec_numbers": attrs.get("ecNumbers", []), "reactants": reactants, "products": products, "raw_reactants": raw_reactants, "raw_products": raw_products, "direction": direction, "enzyme_ids": enzymes, "enzyme_associations": [{"enzyme_id": s.get("subjectId"), "predicate": s.get("predicate"), "sources": s.get("sources", []), "qualifiers": s.get("qualifiers", {})} for s in enzyme_statements], "candidate_proteins": candidate_proteins, "status": "supported" if any((s.get("qualifiers") or {}).get("directExperimentalEvidence") for s in enzyme_statements) else "candidate" if enzymes or candidate_proteins else "unresolved", "carbon_mapping": carbon_mapping, "source_url": entity.get("url"), "directional_rhea_ids": entity.get("identifiers", {}).get("directionalRheaIds", [])})
    for reaction in reactions:
        reaction_class = entities.get(reaction["id"], {}).get("attributes", {}).get("reactionClass", "")
        if str(reaction_class).startswith("non-enzymatic-"):
            reaction["status"] = "non_enzymatic"
    candidate_protein_records = {}
    for item in hypothesis_items:
        for candidate in item.get("candidate_proteins", []):
            p = {**candidate, "proteinId": candidate.get("proteinId") or candidate.get("protein_id"), "geneSymbols": candidate.get("geneSymbols") or candidate.get("gene_symbols", []), "candidateOrigin": candidate.get("candidateOrigin") or candidate.get("candidate_origin")}
            if p.get("proteinId"):
                candidate_protein_records[p["proteinId"]] = p
    for candidates in candidate_by_reaction.values():
        for protein in candidates:
            if protein.get("proteinId"):
                candidate_protein_records.setdefault(protein["proteinId"], enrich_protein(protein))
    for candidates in proteins_by_ec.values():
        for protein in candidates:
            if protein.get("proteinId"):
                candidate_protein_records.setdefault(protein["proteinId"], enrich_protein(protein))
    for protein_id, record in list(candidate_protein_records.items()):
        if protein_id in genome_by_protein:
            candidate_protein_records[protein_id] = {**enrich_protein(record), "sequence_search": genome_by_protein[protein_id].get("sequence_search") or enrich_protein(record).get("sequence_search"), "diamond_search": genome_by_protein[protein_id].get("diamond_search")}
    atom_audit = json.loads(atom_audit_path.read_text()) if atom_audit_path and atom_audit_path.exists() else None
    if atom_audit:
        audit_by_compound = {row["cannabisdb_id"]: row for row in atom_audit.get("compounds", [])}
        for compound in compounds:
            audit_row = audit_by_compound.get(compound.get("id"), {})
            paths = [path for path in audit_row.get("co2_paths", {}).values() if path]
            compound["co2_path_atom_count"] = len(paths)
            compound["co2_path_max_steps"] = max((len(path) for path in paths), default=0)
    lineage_status_counts = {status: sum(target.get("status") == status for target in lineage_targets.values()) for status in ("supported", "candidate", "unresolved")}
    external_id_counts = {db: sum(db in c.get("external_ids", {}) for c in catalog) for db in {db for c in catalog for db in c.get("external_ids", {})}}
    report = {"schema": "cannabis-carbon.networkdb.v1", "sources": {"cannabisdb_compounds": str(compounds_path), "terpedia_network": str(network_path), "identity_crosswalk": str(crosswalk_path), "terpedia_identity_set": str(identity_set_path) if identity_set_path else None, "pubchem_resolution": str(pubchem_path) if pubchem_path else None, "candidate_hypotheses": str(hypotheses_path) if hypotheses_path else None, "genome_search": str(genome_search_path) if genome_search_path else None, "carbon_mapping": str(mapping_path) if mapping_path else None, "carbon_lineage": str(lineage_path) if lineage_path else None, "carbon_atom_audit": str(atom_audit_path) if atom_audit_path else None}, "carbon_source_policy": "CO2 is the only admissible carbon source for Cannabis; no carbon-containing compound is treated as an implicit source.", "compounds": compounds, "reactions": reactions, "identity_links": crosswalk["matches"], "identity_link_candidates": crosswalk.get("candidate_matches", []), "candidate_proteins": list(candidate_protein_records.values()), "candidate_hypotheses": hypothesis_items, "genome_search": genome_search, "carbon_atom_audit": {"source": str(atom_audit_path), "carbon_atoms_total": atom_audit.get("carbon_atoms_total"), "status_counts": atom_audit.get("status_counts")} if atom_audit else None, "coverage": {"cannabisdb_compounds": len(catalog), "terpedia_metabolites": sum(e.get("type") == "metabolite" for e in network["entities"]), "compound_records": len(compounds), "terpedia_reactions": len(reactions), "reaction_records": len(reactions), "exact_identity_links": len(crosswalk["matches"]), "ambiguous_identity_links": crosswalk["ambiguous"], "terpedia_unmatched_metabolites": crosswalk.get("terpedia_unmatched", crosswalk["unmatched"]), "unmatched_cannabisdb_compounds": crosswalk.get("cannabisdb_unmatched", crosswalk["unmatched"]), "connectivity_candidate_identity_links": len(crosswalk.get("candidate_matches", [])), "connectivity_candidate_ambiguous": len(crosswalk.get("candidate_ambiguous_records", [])), "terpedia_identity_set_matches": len(identity_set_by_id), "external_id_coverage": {"records_with_any_external_id": sum(bool(c.get("external_ids")) for c in catalog), "records_without_external_ids": sum(not c.get("external_ids") for c in catalog), "counts_by_database": dict(sorted(external_id_counts.items()))}, "pubchem_resolution": (pubchem or {}).get("summary"), "carbon_mapping_reactions": len(mapping_by_reaction), "carbon_mapping_status_counts": {status: sum(r["carbon_mapping"]["status"] == status for r in reactions) for status in ("inferred", "candidate", "ambiguous", "unresolved", "unavailable")}, "carbon_lineage_edge_totals": lineage_edge_totals, "carbon_lineage_target_status_counts": lineage_status_counts, "candidate_hypotheses": len(hypothesis_items), "candidate_protein_records": len(candidate_protein_records), "candidate_proteins_with_sequence": sum(bool(p.get("sequence_search", {}).get("sequence_present")) for p in candidate_protein_records.values()), "candidate_proteins_with_diamond_hits": sum(bool((p.get("diamond_search") or {}).get("hits")) for p in candidate_protein_records.values()), "candidate_proteins_meeting_strong_hit_threshold": sum(bool((p.get("diamond_search") or {}).get("strong_candidate_hit")) for p in candidate_protein_records.values()), "reactions_with_candidate_proteins": sum(bool(r["candidate_proteins"]) for r in reactions), "reactions_with_exact_ec_annotation_candidates": sum(any(p.get("candidateOrigin") == "Terpedia reaction EC to protein exact-EC join" for p in r["candidate_proteins"]) for r in reactions)}, "claim_boundary": "This is a unified inventory and reaction database. Presence of a compound, reaction, enzyme association, or identity link does not establish in-vivo cannabis biosynthesis; candidate proteins remain hypotheses and sequence similarity is not functional proof."}
    output.parent.mkdir(parents=True, exist_ok=True)
    report["hypothetical_connections"] = hypothetical_connections
    report["sources"]["terpedia_hypothetical_forward_connections"] = str(hypothetical_connections_path) if hypothetical_connections_path else None
    report["coverage"]["hypothetical_forward_connections"] = len(hypothetical_connections)
    report["coverage"]["hypothetical_identity_compounds"] = sum(compound.get("namespace") == "terpedia_identity_set" for compound in compounds)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return report["coverage"]
