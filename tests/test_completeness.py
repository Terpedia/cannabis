import json
import gzip

from cannabis_carbon.completeness import compute_completeness


def test_completeness_separates_missing_metabolites_and_enzymes(tmp_path):
    network = {"entities": [
        {"id": "m:a", "type": "metabolite"}, {"id": "m:b", "type": "metabolite"},
        {"id": "m:c", "type": "metabolite"}, {"id": "r:1", "type": "biochemical_reaction"}
    ], "statements": [
        {"subjectId": "r:1", "predicate": "has_reactant", "objectEntityId": "m:a"},
        {"subjectId": "r:1", "predicate": "has_product", "objectEntityId": "m:b"}
    ]}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": [{"carbon_atom_count": 1}]}))
    result = compute_completeness(network_path, compounds_path)
    assert result["terpedia"]["metabolites_without_reactions"] == 1
    assert result["terpedia"]["reactions_without_enzyme_association"] == 1


def test_has_catalytic_activity_counts_as_enzyme_association(tmp_path):
    network = {"entities": [{"id": "r:1", "type": "biochemical_reaction"}], "statements": [
        {"subjectId": "p:1", "predicate": "has_catalytic_activity", "objectEntityId": "r:1"}
    ]}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    result = compute_completeness(network_path, compounds_path)
    assert result["terpedia"]["reactions_with_enzyme_association"] == 1
    assert result["terpedia"]["reactions_without_enzyme_association"] == 0


def test_non_enzymatic_reactions_are_not_enzyme_gaps(tmp_path):
    network = {"entities": [{"id": "r:1", "type": "biochemical_reaction", "attributes": {"reactionClass": "non-enzymatic-decarboxylation"}}], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    result = compute_completeness(network_path, compounds_path)
    assert result["terpedia"]["non_enzymatic_reactions"] == 1
    assert result["terpedia"]["enzyme_requiring_reactions"] == 0
    assert result["terpedia"]["reactions_without_enzyme_association"] == 0


def test_completeness_can_include_co2_lineage(tmp_path):
    network = {"entities": [], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    lineage_path = tmp_path / "lineage.json"
    lineage_path.write_text(json.dumps({"target_summary": {"supported": 1}, "reachable_carbon_nodes": 2, "resolved_carbon_edges": 1, "inferred_carbon_edges": 1, "candidate_carbon_edges": 0, "external_carbon_input_entity_count": 3, "carbon_source_policy": "CO2 only"}))
    result = compute_completeness(network_path, compounds_path, lineage_path=lineage_path)
    assert result["coverage"]["co2_lineage"]["target_summary"]["supported"] == 1


def test_completeness_reports_mapping_blockers(tmp_path):
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump({"entities": [], "statements": []}, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    mapping_path = tmp_path / "mapping.json"
    mapping_path.write_text(json.dumps({"carbon_counts": {"mapped_carbon_atoms": 1, "unresolved_or_ambiguous_carbon_atoms": 2, "product_carbon_atoms": 3, "mapping_coverage_percent": 33.3}, "status_counts": {}, "reactions": [{"reaction_id": "r1", "mappings": [{"status": "ambiguous"}, {"status": "unresolved"}, {"status": "inferred"}]}]}))
    result = compute_completeness(network_path, compounds_path, mapping_path=mapping_path)
    assert result["coverage"]["carbon_mapping_blockers"] == {"reactions_with_blocked_product_carbon_rows": 1, "product_carbon_row_status_counts": {"ambiguous": 1, "inferred": 1, "unresolved": 1}, "blocked_product_carbon_rows": 2}


def test_completeness_reports_atom_percentages(tmp_path):
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump({"entities": [], "statements": []}, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    audit_path = tmp_path / "audit.json"
    audit_path.write_text(json.dumps({"carbon_atoms_total": 10, "compound_count": 1, "status_counts": {"supported": 1, "candidate": 2, "inferred": 3, "unresolved": 4}}))
    result = compute_completeness(network_path, compounds_path, atom_audit_path=audit_path)
    assert result["coverage"]["carbon_atom_audit"]["evidence_bearing_percent"] == 60.0
    assert result["coverage"]["carbon_atom_audit"]["unresolved_percent"] == 40.0


def test_completeness_distinguishes_nonexact_from_unresolved_identity(tmp_path):
    network = {"entities": [], "statements": []}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": [{"id": "CDB1", "carbon_atom_count": 1}, {"id": "CDB2", "carbon_atom_count": 1}]}))
    crosswalk = tmp_path / "crosswalk.json"
    crosswalk.write_text(json.dumps({"matches": [{"cannabisdb": {"cannabisdb_id": "CDB1"}}], "ambiguous": 0, "unmatched": 1, "cannabisdb_unmatched": 1, "connectivity_candidate_matches": 1, "connectivity_candidate_ambiguous": 0, "tautomer_candidate_matches": 1, "tautomer_candidate_ambiguous": 0}))
    result = compute_completeness(network_path, compounds_path, crosswalk_path=crosswalk)
    assert result["cannabisdb"]["compounds_without_exact_terpedia_identity"] == 1
    assert result["cannabisdb"]["compounds_without_any_identity_resolution"] == 1
    assert result["cannabisdb"]["tautomer_candidate_identity_links"] == 1


def test_completeness_splits_unannotated_reactions_by_candidate_proteins(tmp_path):
    network = {"entities": [
        {"id": "m:a", "type": "metabolite"}, {"id": "m:b", "type": "metabolite"},
        {"id": "r:1", "type": "biochemical_reaction"},
        {"id": "r:2", "type": "biochemical_reaction"},
    ], "statements": [
        {"subjectId": "r:1", "predicate": "has_reactant", "objectEntityId": "m:a"},
        {"subjectId": "r:1", "predicate": "has_product", "objectEntityId": "m:b"},
        {"subjectId": "r:2", "predicate": "has_reactant", "objectEntityId": "m:a"},
        {"subjectId": "r:2", "predicate": "has_product", "objectEntityId": "m:b"},
    ]}
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump(network, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    hypotheses_path = tmp_path / "hypotheses.json"
    hypotheses_path.write_text(json.dumps({"items": [{"reaction_id": "r:1", "candidate_proteins": [{"id": "p:1"}]}]}))
    result = compute_completeness(network_path, compounds_path, hypotheses_path=hypotheses_path)
    assert result["terpedia"]["reactions_without_enzyme_association"] == 2
    assert result["terpedia"]["reactions_without_enzyme_with_candidate_proteins"] == 1
    assert result["terpedia"]["reactions_without_enzyme_without_candidate_proteins"] == 1


def test_completeness_reports_co2_target_triage(tmp_path):
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump({"entities": [], "statements": []}, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": [{"id": "CDB1", "label": "Cannabifoo", "aliases": [], "carbon_atom_count": 3}]}))
    lineage_path = tmp_path / "lineage.json"
    lineage_path.write_text(json.dumps({"target_summary": {"supported": 0, "candidate": 1, "unresolved": 0}, "reachable_carbon_nodes": 1, "resolved_carbon_edges": 1, "inferred_carbon_edges": 1, "candidate_carbon_edges": 0, "external_carbon_input_entity_count": 0, "carbon_source_policy": "CO2 only", "targets": [{"cannabisdb_id": "CDB1", "status": "candidate", "identity_status": "exact", "carbon_atom_count": 3, "reason": "partial"}]}))
    result = compute_completeness(network_path, compounds_path, lineage_path=lineage_path)
    triage = result["coverage"]["co2_lineage"]["target_triage"]
    assert triage["target_counts_by_status_and_identity"] == {"candidate:exact": 1}
    assert triage["carbon_atoms_by_target_status"] == {"candidate": 3}
    assert triage["specialty_target_counts_by_status"] == {"candidate": 1}


def test_completeness_reports_hypothesis_layer_separately(tmp_path):
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump({"entities": [], "statements": []}, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    networkdb_path = tmp_path / "networkdb.json"
    networkdb_path.write_text(json.dumps({"compounds": [{"id": "x"}], "hypothetical_connections": [{"status": "candidate", "enzyme_evidence": [{}]}, {"status": "unresolved"}], "coverage": {"hypothetical_reaction_inventory": 3, "hypothetical_product_inventory": 2, "hypothetical_missing_substrate_nodes": 1}}))
    result = compute_completeness(network_path, compounds_path, networkdb_path=networkdb_path)
    assert result["hypothesis_layer"]["hypothesis_edges"] == 2
    assert result["hypothesis_layer"]["hypothesis_edges_with_enzyme_evidence"] == 1
    assert result["hypothesis_layer"]["hypothesis_reaction_inventory"] == 3


def test_completeness_reports_candidate_path_carbon_layer(tmp_path):
    network_path = tmp_path / "network.json.gz"
    with gzip.open(network_path, "wt") as handle: json.dump({"entities": [], "statements": []}, handle)
    compounds_path = tmp_path / "compounds.json"
    compounds_path.write_text(json.dumps({"compounds": []}))
    candidate_path = tmp_path / "candidate-path-carbon.json"
    candidate_path.write_text(json.dumps({
        "path_count": 2,
        "paths_with_complete_product_carbon_mapping": 1,
        "mapped_product_carbon_atoms": 4,
        "unresolved_product_carbon_atoms": 1,
        "rows": [
            {"candidate_product_terpene_id": "T1", "carbon_mapping": {"status": "inferred"}, "core_path_carbon_edges": [{"from_atom": None}]},
            {"candidate_product_terpene_id": "T2", "carbon_mapping": {"status": "unresolved"}, "core_path_carbon_edges": []},
        ],
    }))
    result = compute_completeness(network_path, compounds_path, candidate_path_carbon_path=candidate_path)
    layer = result["coverage"]["candidate_co2_path_carbon_layer"]
    assert layer["path_count"] == 2
    assert layer["candidate_product_count"] == 2
    assert layer["carbon_mapping_status_counts"] == {"inferred": 1, "unresolved": 1}
    assert layer["atom_level_core_path_edges"] == 1
