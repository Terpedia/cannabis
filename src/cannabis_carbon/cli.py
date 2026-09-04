from __future__ import annotations

import argparse
import hashlib
import json
import ssl
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from .candidates import CandidateEvidence, rank_candidate
from .ingest import ingest_sdf
from .terpedia import cytoscape_elements, load_network
from .reaction_report import build_reaction_report
from .completeness import compute_completeness
from .hypotheses import build_candidate_queue, build_carbon_mapping_queue
from .crosswalk import build_crosswalk
from .balance import audit_balances
from .hypothesis_balance import audit_hypothesis_balances
from .hypothesis_mapping import build_hypothesis_mapping
from .hypothesis_lineage import build_hypothesis_lineage
from .lineage import build_carbon_lineage
from .networkdb import build_networkdb, build_map_snapshot
from .genome import build_genome_search
from .inventory import build_specialty_inventory
from .test_hypotheses import build_test_hypotheses
from .validate import validate_artifacts
from .pubchem import resolve_pubchem
from .cannabisdb_xml import enrich_compounds_with_xrefs, extract_terpedia_table
from .identity_set import refresh_identity_set
from .identity_set_paths import refresh_identity_set_upstream, refresh_identity_set_connectivity_upstream, refresh_identity_set_candidate_expansion, build_candidate_expansion_bridges, build_candidate_expansion_carbon_mapping, map_identity_set_upstream, build_identity_set_core_bridges
from .pubchem_xrefs import retrieve_pubchem_chebi_xrefs

DOWNLOADS = {
    "compounds.sdf": "https://cannabisdatabase.ca/simple/download_compound_as_sdf",
    "compounds.xml": "https://cannabisdatabase.ca/simple/download_compound_as_xml",
    "proteins.xml": "https://cannabisdatabase.ca/simple/download_protein_as_xml",
}


def download(out: Path, insecure: bool = False) -> None:
    out.mkdir(parents=True, exist_ok=True)
    context = ssl._create_unverified_context() if insecure else ssl.create_default_context()
    manifest = {"source": "Cannabis Compound Database", "retrieved_at": datetime.now(timezone.utc).isoformat(), "files": []}
    for name, url in DOWNLOADS.items():
        target = out / name
        with urllib.request.urlopen(url, context=context) as response, target.open("wb") as handle:
            while chunk := response.read(1024 * 1024):
                handle.write(chunk)
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        manifest["files"].append({"name": name, "url": url, "sha256": digest, "bytes": target.stat().st_size})
        print(f"{name}: {target.stat().st_size:,} bytes sha256={digest}")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def inspect_sdf(path: Path) -> None:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    count = valid = carbon = 0
    for mol in supplier:
        count += 1
        if mol is None:
            continue
        valid += 1
        carbon += sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    print(json.dumps({"file": str(path), "records": count, "valid_structures": valid, "carbon_atoms": carbon}, indent=2))


def rank_candidates(path: Path) -> None:
    record = json.loads(path.read_text())
    candidate = CandidateEvidence(
        protein_id=record["protein_id"], reaction_id=record["reaction_id"],
        identity=record.get("identity"), coverage=record.get("coverage"),
        profile_score=record.get("profile_score"),
        catalytic_motif=record.get("catalytic_motif", False),
        complete_domains=record.get("complete_domains", False),
        localization_support=record.get("localization_support", False),
        expression_support=record.get("expression_support", False),
        source_urls=tuple(record.get("source_urls", [])),
    )
    print(json.dumps(rank_candidate(candidate), indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p_download = sub.add_parser("download")
    p_download.add_argument("--out", type=Path, default=Path("data/raw"))
    p_download.add_argument("--insecure-download", action="store_true", help="Allow the current CannabisDB expired TLS certificate")
    p_inspect = sub.add_parser("inspect-sdf")
    p_inspect.add_argument("path", type=Path)
    p_rank = sub.add_parser("rank-candidate")
    p_rank.add_argument("path", type=Path)
    p_ingest = sub.add_parser("ingest-sdf")
    p_ingest.add_argument("source", type=Path)
    p_ingest.add_argument("--graph-out", type=Path, default=Path("docs/data/compounds.json"))
    p_ingest.add_argument("--report-out", type=Path, default=Path("data/reports/carbon-coverage.json"))
    p_export = sub.add_parser("export-terpedia-graph")
    p_export.add_argument("source", type=Path)
    p_export.add_argument("--out", type=Path, default=Path("docs/data/terpedia-network.json"))
    p_map = sub.add_parser("map-reactions")
    p_map.add_argument("source", type=Path)
    p_map.add_argument("--out", type=Path, default=Path("data/reports/reaction-carbon-mapping.json"))
    p_complete = sub.add_parser("completeness")
    p_complete.add_argument("network", type=Path)
    p_complete.add_argument("compounds", type=Path, default=Path("docs/data/compounds.json"), nargs="?")
    p_complete.add_argument("--mapping", type=Path)
    p_complete.add_argument("--crosswalk", type=Path)
    p_complete.add_argument("--lineage", type=Path)
    p_complete.add_argument("--atom-audit", type=Path, default=Path("data/reports/carbon-atom-audit.json"))
    p_complete.add_argument("--hypotheses", type=Path, default=Path("data/reports/candidate-work-queue.json"))
    p_complete.add_argument("--pubchem", type=Path, default=Path("data/reports/pubchem-resolution.json"))
    p_complete.add_argument("--networkdb", type=Path, default=Path("docs/data/networkdb.json"))
    p_complete.add_argument("--hypothesis-lineage", type=Path, default=Path("data/reports/hypothesis-lineage.json"))
    p_complete.add_argument("--candidate-path-carbon", type=Path, default=Path("data/reports/terpene-identity-set-reversible-candidate-lineage-carbon.json"))
    p_complete.add_argument("--out", type=Path, default=Path("data/reports/completeness.json"))
    p_queue = sub.add_parser("candidate-queue")
    p_queue.add_argument("source", type=Path)
    p_queue.add_argument("--out", type=Path, default=Path("data/reports/candidate-work-queue.json"))
    p_mapping_queue = sub.add_parser("carbon-mapping-queue")
    p_mapping_queue.add_argument("mapping", type=Path)
    p_mapping_queue.add_argument("networkdb", type=Path)
    p_mapping_queue.add_argument("--out", type=Path, default=Path("data/reports/carbon-mapping-work-queue.json"))
    p_crosswalk = sub.add_parser("crosswalk")
    p_crosswalk.add_argument("cannabisdb_sdf", type=Path)
    p_crosswalk.add_argument("terpedia_network", type=Path)
    p_crosswalk.add_argument("--compounds", type=Path, help="Normalized CannabisDB catalog with XML external identifiers")
    p_crosswalk.add_argument("--pubchem-chebi", type=Path, help="PubChem ChEBI cross-reference report")
    p_crosswalk.add_argument("--out", type=Path, default=Path("data/reports/identity-crosswalk.json"))
    p_balance = sub.add_parser("balance-audit")
    p_balance.add_argument("network", type=Path)
    p_balance.add_argument("--out", type=Path, default=Path("data/reports/phase1-balance-audit.json"))
    p_hypothesis_balance = sub.add_parser("hypothesis-balance-audit")
    p_hypothesis_balance.add_argument("source", type=Path, default=Path("data/terpedia/hypothetical-forward-connections.json"), nargs="?")
    p_hypothesis_balance.add_argument("--out", type=Path, default=Path("data/reports/terpedia-hypothesis-balance-audit.json"))
    p_candidate_balance = sub.add_parser("candidate-expansion-balance-audit")
    p_candidate_balance.add_argument("source", type=Path, default=Path("data/reports/terpene-identity-set-candidate-expansion.json"), nargs="?")
    p_candidate_balance.add_argument("--out", type=Path, default=Path("data/reports/terpene-identity-set-candidate-expansion-balance-audit.json"))
    p_hypothesis_mapping = sub.add_parser("hypothesis-carbon-mapping")
    p_hypothesis_mapping.add_argument("source", type=Path, default=Path("data/terpedia/hypothetical-forward-connections.json"), nargs="?")
    p_hypothesis_mapping.add_argument("--out", type=Path, default=Path("data/reports/terpedia-hypothesis-carbon-mapping.json"))
    p_lineage = sub.add_parser("carbon-lineage")
    p_lineage.add_argument("network", type=Path)
    p_lineage.add_argument("mapping", type=Path)
    p_lineage.add_argument("crosswalk", type=Path)
    p_lineage.add_argument("compounds", type=Path)
    p_lineage.add_argument("--directions", type=Path, default=Path("data/terpedia/directional-reaction-overrides.json"))
    p_lineage.add_argument("--out", type=Path, default=Path("data/reports/carbon-lineage.json"))
    p_atom_audit = sub.add_parser("carbon-atom-audit")
    p_atom_audit.add_argument("network", type=Path)
    p_atom_audit.add_argument("lineage", type=Path)
    p_atom_audit.add_argument("crosswalk", type=Path)
    p_atom_audit.add_argument("compounds", type=Path)
    p_atom_audit.add_argument("--networkdb", type=Path, default=Path("data/reports/networkdb.json"))
    p_atom_audit.add_argument("--out", type=Path, default=Path("data/reports/carbon-atom-audit.json"))
    p_networkdb = sub.add_parser("networkdb")
    p_networkdb.add_argument("network", type=Path)
    p_networkdb.add_argument("compounds", type=Path)
    p_networkdb.add_argument("crosswalk", type=Path)
    p_networkdb.add_argument("--hypotheses", type=Path, default=Path("data/reports/candidate-work-queue.json"))
    p_networkdb.add_argument("--genome-search", type=Path, default=Path("data/reports/genome-candidate-search.json"))
    p_networkdb.add_argument("--genome-fasta", type=Path, default=Path("data/raw/UP000583929.fasta"))
    p_networkdb.add_argument("--mapping", type=Path, default=Path("data/reports/reaction-carbon-mapping.json"))
    p_networkdb.add_argument("--lineage", type=Path, default=Path("data/reports/carbon-lineage.json"))
    p_networkdb.add_argument("--atom-audit", type=Path, default=Path("data/reports/carbon-atom-audit.json"))
    p_networkdb.add_argument("--pubchem", type=Path, default=Path("data/reports/pubchem-resolution.json"))
    p_networkdb.add_argument("--identity-set", type=Path, default=Path("data/reports/terpene-identity-set-match.json"))
    p_networkdb.add_argument("--hypothetical-connections", type=Path, default=Path("data/terpedia/hypothetical-forward-connections.json"))
    p_networkdb.add_argument("--hypothetical-reactions", type=Path, default=Path("data/terpedia/hypothetical-reaction-inventory.json"))
    p_networkdb.add_argument("--hypothesis-enzyme-evidence", type=Path, default=Path("data/terpedia/terpene-enzyme-reaction-gene-evidence.json"))
    p_networkdb.add_argument("--hypothesis-enzyme-catalog", type=Path, default=Path("data/terpedia/terpene-biotransformation-enzyme-catalog.json"))
    p_networkdb.add_argument("--hypothesis-balance", type=Path, default=Path("data/reports/terpedia-hypothesis-balance-audit.json"))
    p_networkdb.add_argument("--identity-set-bridges", type=Path, default=Path("data/reports/terpedia-identity-set-core-bridges.json"))
    p_networkdb.add_argument("--identity-set-connectivity-upstream", type=Path, default=Path("data/reports/terpene-identity-set-connectivity-upstream.json"))
    p_networkdb.add_argument("--out", type=Path, default=Path("docs/data/networkdb.json"))
    p_map_snapshot = sub.add_parser("map-snapshot")
    p_map_snapshot.add_argument("networkdb", type=Path)
    p_map_snapshot.add_argument("--lineage", type=Path, default=Path("data/reports/carbon-lineage.json"))
    p_map_snapshot.add_argument("--focus-out", type=Path)
    p_map_snapshot.add_argument("--out", type=Path, default=Path("docs/data/network-map.json"))
    p_enzyme_gaps = sub.add_parser("enzyme-gap-audit")
    p_enzyme_gaps.add_argument("networkdb", type=Path, default=Path("docs/data/networkdb.json"), nargs="?")
    p_enzyme_gaps.add_argument("--out", type=Path, default=Path("data/reports/enzyme-gap-audit.json"))
    p_hyp_lineage = sub.add_parser("hypothesis-lineage")
    p_hyp_lineage.add_argument("networkdb", type=Path, default=Path("docs/data/networkdb.json"), nargs="?")
    p_hyp_lineage.add_argument("--out", type=Path, default=Path("data/reports/hypothesis-lineage.json"))
    p_genome = sub.add_parser("genome-search")
    p_genome.add_argument("queue", type=Path)
    p_genome.add_argument("fasta", type=Path)
    p_genome.add_argument("--diamond-hits", type=Path)
    p_genome.add_argument("--reference-tsv", type=Path)
    p_genome.add_argument("--out", type=Path, default=Path("data/reports/genome-candidate-search.json"))
    p_inventory = sub.add_parser("specialty-inventory")
    p_inventory.add_argument("compounds", type=Path)
    p_inventory.add_argument("crosswalk", type=Path)
    p_inventory.add_argument("network", type=Path)
    p_inventory.add_argument("--lineage", type=Path, default=Path("data/reports/carbon-lineage.json"))
    p_inventory.add_argument("--pubchem", type=Path, default=Path("data/reports/pubchem-resolution.json"))
    p_inventory.add_argument("--out", type=Path, default=Path("data/reports/named-specialty-inventory.json"))
    p_tests = sub.add_parser("test-hypotheses")
    p_tests.add_argument("queue", type=Path)
    p_tests.add_argument("network", type=Path)
    p_tests.add_argument("--lineage", type=Path, default=Path("data/reports/carbon-lineage.json"))
    p_tests.add_argument("--compounds", type=Path, default=Path("docs/data/compounds.json"))
    p_tests.add_argument("--out", type=Path, default=Path("data/reports/testable-hypotheses.json"))
    p_validate = sub.add_parser("validate-artifacts")
    p_validate.add_argument("atom_audit", type=Path)
    p_validate.add_argument("mapping", type=Path)
    p_validate.add_argument("balance", type=Path)
    p_validate.add_argument("compounds", type=Path)
    p_validate.add_argument("--mapping-queue", type=Path, default=Path("data/reports/carbon-mapping-work-queue.json"))
    p_validate.add_argument("--out", type=Path, default=Path("data/reports/artifact-validation.json"))
    p_pubchem = sub.add_parser("pubchem-resolve")
    p_pubchem.add_argument("compounds", type=Path, default=Path("docs/data/compounds.json"), nargs="?")
    p_pubchem.add_argument("--out", type=Path, default=Path("data/reports/pubchem-resolution.json"))
    p_pubchem.add_argument("--batch-size", type=int, default=25)
    p_pubchem.add_argument("--pause", type=float, default=0.25)
    p_pubchem.add_argument("--workers", type=int, default=4)
    p_pubchem.add_argument("--cache", type=Path, default=Path("data/reports/pubchem-cache.json"))
    p_pubchem.add_argument("--method", choices=("bulk", "batch"), default="bulk")
    p_pubchem_xrefs = sub.add_parser("pubchem-chebi-xrefs")
    p_pubchem_xrefs.add_argument("pubchem", type=Path, default=Path("data/reports/pubchem-resolution.json"), nargs="?")
    p_pubchem_xrefs.add_argument("--out", type=Path, default=Path("data/reports/pubchem-chebi-xrefs.json"))
    p_pubchem_xrefs.add_argument("--workers", type=int, default=4)
    p_pubchem_xrefs.add_argument("--pause", type=float, default=0.15)
    p_xrefs = sub.add_parser("enrich-cannabisdb-xrefs")
    p_xrefs.add_argument("xml", type=Path)
    p_xrefs.add_argument("compounds", type=Path, default=Path("docs/data/compounds.json"), nargs="?")
    p_xrefs.add_argument("--out", type=Path, default=Path("docs/data/compounds.json"))
    p_xrefs.add_argument("--report", type=Path, default=Path("data/reports/cannabisdb-xrefs.json"))
    p_identity_refresh = sub.add_parser("refresh-identity-set")
    p_identity_refresh.add_argument("compounds", type=Path, default=Path("docs/data/compounds.json"), nargs="?")
    p_identity_refresh.add_argument("--out", type=Path, default=Path("data/reports/terpene-identity-set-match.json"))
    p_identity_refresh.add_argument("--bq", default="bq")
    p_identity_upstream = sub.add_parser("refresh-identity-set-upstream")
    p_identity_upstream.add_argument("compounds", type=Path, default=Path("docs/data/compounds.json"), nargs="?")
    p_identity_upstream.add_argument("--out", type=Path, default=Path("data/reports/terpedia-identity-set-upstream.json"))
    p_identity_upstream.add_argument("--bq", default="bq")
    p_identity_connectivity_upstream = sub.add_parser("refresh-identity-set-connectivity-upstream")
    p_identity_connectivity_upstream.add_argument("compounds", type=Path, default=Path("docs/data/compounds.json"), nargs="?")
    p_identity_connectivity_upstream.add_argument("identity_set", type=Path, default=Path("data/reports/terpene-identity-set-match.json"), nargs="?")
    p_identity_connectivity_upstream.add_argument("--out", type=Path, default=Path("data/reports/terpene-identity-set-connectivity-upstream.json"))
    p_identity_connectivity_upstream.add_argument("--bq", default="bq")
    p_identity_candidate_expansion = sub.add_parser("refresh-identity-set-candidate-expansion")
    p_identity_candidate_expansion.add_argument("connectivity", type=Path, default=Path("data/reports/terpene-identity-set-connectivity-upstream.json"), nargs="?")
    p_identity_candidate_expansion.add_argument("--out", type=Path, default=Path("data/reports/terpene-identity-set-candidate-expansion.json"))
    p_identity_candidate_expansion.add_argument("--depth", type=int, default=3)
    p_identity_candidate_expansion.add_argument("--bq", default="bq")
    p_identity_candidate_bridges = sub.add_parser("candidate-expansion-bridges")
    p_identity_candidate_bridges.add_argument("expansion", type=Path, default=Path("data/reports/terpene-identity-set-candidate-expansion.json"), nargs="?")
    p_identity_candidate_bridges.add_argument("network", type=Path, default=Path("data/terpedia/cannabis-sativa-metabolic-network.json.gz"), nargs="?")
    p_identity_candidate_bridges.add_argument("--lineage", type=Path, default=Path("data/reports/carbon-lineage.json"))
    p_identity_candidate_bridges.add_argument("--out", type=Path, default=Path("data/reports/terpene-identity-set-candidate-expansion-bridges.json"))
    p_identity_candidate_mapping = sub.add_parser("candidate-expansion-carbon-mapping")
    p_identity_candidate_mapping.add_argument("expansion", type=Path, default=Path("data/reports/terpene-identity-set-candidate-expansion.json"), nargs="?")
    p_identity_candidate_mapping.add_argument("bridges", type=Path, default=Path("data/reports/terpene-identity-set-candidate-expansion-bridges.json"), nargs="?")
    p_identity_candidate_mapping.add_argument("--out", type=Path, default=Path("data/reports/terpene-identity-set-candidate-expansion-carbon-mapping.json"))
    p_candidate_lineage = sub.add_parser("reversible-candidate-lineage")
    p_candidate_lineage.add_argument("bridges", type=Path, default=Path("data/reports/terpene-identity-set-candidate-expansion-bridges.json"), nargs="?")
    p_candidate_lineage.add_argument("lineage", type=Path, default=Path("data/reports/carbon-lineage.json"), nargs="?")
    p_candidate_lineage.add_argument("--out", type=Path, default=Path("data/reports/terpene-identity-set-reversible-candidate-lineage.json"))
    p_candidate_lineage_carbon = sub.add_parser("reversible-candidate-lineage-carbon")
    p_candidate_lineage_carbon.add_argument("lineage", type=Path, default=Path("data/reports/terpene-identity-set-reversible-candidate-lineage.json"), nargs="?")
    p_candidate_lineage_carbon.add_argument("mapping", type=Path, default=Path("data/reports/terpene-identity-set-candidate-expansion-carbon-mapping.json"), nargs="?")
    p_candidate_lineage_carbon.add_argument("--out", type=Path, default=Path("data/reports/terpene-identity-set-reversible-candidate-lineage-carbon.json"))
    p_identity_map = sub.add_parser("map-identity-set-upstream")
    p_identity_map.add_argument("source", type=Path, default=Path("data/reports/terpedia-identity-set-upstream.json"), nargs="?")
    p_identity_map.add_argument("--out", type=Path, default=Path("data/reports/terpedia-identity-set-upstream-mapped.json"))
    p_identity_bridges = sub.add_parser("identity-set-core-bridges")
    p_identity_bridges.add_argument("source", type=Path, default=Path("data/reports/terpedia-identity-set-upstream-mapped.json"), nargs="?")
    p_identity_bridges.add_argument("network", type=Path, default=Path("data/terpedia/cannabis-sativa-metabolic-network.json.gz"), nargs="?")
    p_identity_bridges.add_argument("--out", type=Path, default=Path("data/reports/terpedia-identity-set-core-bridges.json"))
    p_table = sub.add_parser("extract-cannabisdb-table")
    p_table.add_argument("xml", type=Path)
    p_table.add_argument("--out", type=Path, default=Path("data/terpedia/cannabisdb-compounds.json"))
    p_table.add_argument("--report", type=Path, default=Path("data/reports/cannabisdb-table.json"))
    args = parser.parse_args()
    if args.command == "download":
        download(args.out, args.insecure_download)
    elif args.command == "inspect-sdf":
        inspect_sdf(args.path)
    elif args.command == "rank-candidate":
        rank_candidates(args.path)
    elif args.command == "ingest-sdf":
        print(json.dumps(ingest_sdf(args.source, args.graph_out, args.report_out), indent=2))
    elif args.command == "export-terpedia-graph":
        directions_path = args.source.parent / "directional-reaction-overrides.json"
        directions = json.loads(directions_path.read_text()) if directions_path.exists() else {}
        graph = cytoscape_elements(load_network(args.source), directions)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(graph, separators=(",", ":")) + "\n")
        print(json.dumps(graph["stats"], indent=2))
    elif args.command == "map-reactions":
        print(json.dumps(build_reaction_report(args.source, args.out), indent=2))
    elif args.command == "completeness":
        result = compute_completeness(args.network, args.compounds, args.mapping, args.crosswalk, args.lineage, args.atom_audit, args.hypotheses, args.pubchem, args.networkdb, args.hypothesis_lineage, args.candidate_path_carbon)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
    elif args.command == "candidate-queue":
        print(json.dumps(build_candidate_queue(args.source, args.out), indent=2))
    elif args.command == "carbon-mapping-queue":
        print(json.dumps(build_carbon_mapping_queue(args.mapping, args.networkdb, args.out), indent=2))
    elif args.command == "crosswalk":
        print(json.dumps(build_crosswalk(args.cannabisdb_sdf, args.terpedia_network, args.out, args.compounds, args.pubchem_chebi), indent=2))
    elif args.command == "balance-audit":
        print(json.dumps(audit_balances(args.network, args.out), indent=2))
    elif args.command == "hypothesis-balance-audit":
        print(json.dumps(audit_hypothesis_balances(args.source, args.out), indent=2))
    elif args.command == "candidate-expansion-balance-audit":
        from .candidate_balance import audit_candidate_expansion_balances
        print(json.dumps(audit_candidate_expansion_balances(args.source, args.out), indent=2))
    elif args.command == "hypothesis-carbon-mapping":
        print(json.dumps(build_hypothesis_mapping(args.source, args.out), indent=2))
    elif args.command == "carbon-lineage":
        print(json.dumps(build_carbon_lineage(args.network, args.mapping, args.crosswalk, args.compounds, args.out, args.directions), indent=2))
    elif args.command == "carbon-atom-audit":
        from .lineage import build_carbon_atom_audit
        print(json.dumps(build_carbon_atom_audit(args.network, args.lineage, args.crosswalk, args.compounds, args.out, args.networkdb), indent=2))
    elif args.command == "networkdb":
        print(json.dumps(build_networkdb(args.network, args.compounds, args.crosswalk, args.out, args.hypotheses, args.genome_search, args.genome_fasta, args.mapping, args.lineage, args.atom_audit, args.pubchem, args.identity_set, args.hypothetical_connections, args.hypothetical_reactions, args.hypothesis_enzyme_evidence, args.hypothesis_enzyme_catalog, args.hypothesis_balance, args.identity_set_bridges, args.identity_set_connectivity_upstream), indent=2))
    elif args.command == "map-snapshot":
        print(json.dumps(build_map_snapshot(args.networkdb, args.out, args.lineage, args.focus_out), indent=2))
    elif args.command == "enzyme-gap-audit":
        from .enzyme_gaps import build_enzyme_gap_audit
        print(json.dumps(build_enzyme_gap_audit(args.networkdb, args.out), indent=2))
    elif args.command == "hypothesis-lineage":
        print(json.dumps(build_hypothesis_lineage(args.networkdb, args.out), indent=2))
    elif args.command == "genome-search":
        print(json.dumps(build_genome_search(args.queue, args.fasta, args.out, args.diamond_hits, args.reference_tsv), indent=2))
    elif args.command == "specialty-inventory":
        print(json.dumps(build_specialty_inventory(args.compounds, args.crosswalk, args.network, args.out, args.lineage, args.pubchem), indent=2))
    elif args.command == "test-hypotheses":
        print(json.dumps(build_test_hypotheses(args.queue, args.network, args.out, args.lineage, args.compounds), indent=2))
    elif args.command == "validate-artifacts":
        result = validate_artifacts(args.atom_audit, args.mapping, args.balance, args.compounds, args.out, args.mapping_queue)
        print(json.dumps(result, indent=2))
        if not result["valid"]:
            raise SystemExit(1)
    elif args.command == "pubchem-resolve":
        print(json.dumps(resolve_pubchem(args.compounds, args.out, args.batch_size, args.pause, args.workers, args.cache, args.method), indent=2))
    elif args.command == "pubchem-chebi-xrefs":
        print(json.dumps(retrieve_pubchem_chebi_xrefs(args.pubchem, args.out, args.workers, args.pause), indent=2))
    elif args.command == "enrich-cannabisdb-xrefs":
        print(json.dumps(enrich_compounds_with_xrefs(args.xml, args.compounds, args.out, args.report), indent=2))
    elif args.command == "refresh-identity-set":
        print(json.dumps(refresh_identity_set(args.compounds, args.out, args.bq), indent=2))
    elif args.command == "refresh-identity-set-upstream":
        print(json.dumps(refresh_identity_set_upstream(args.compounds, args.out, args.bq), indent=2))
    elif args.command == "refresh-identity-set-connectivity-upstream":
        print(json.dumps(refresh_identity_set_connectivity_upstream(args.compounds, args.identity_set, args.out, args.bq), indent=2))
    elif args.command == "refresh-identity-set-candidate-expansion":
        print(json.dumps(refresh_identity_set_candidate_expansion(args.connectivity, args.out, args.bq, args.depth), indent=2))
    elif args.command == "candidate-expansion-bridges":
        print(json.dumps(build_candidate_expansion_bridges(args.expansion, args.network, args.lineage, args.out), indent=2))
    elif args.command == "candidate-expansion-carbon-mapping":
        print(json.dumps(build_candidate_expansion_carbon_mapping(args.expansion, args.bridges, args.out), indent=2))
    elif args.command == "reversible-candidate-lineage":
        from .candidate_lineage import build_reversible_candidate_lineage
        print(json.dumps(build_reversible_candidate_lineage(args.bridges, args.lineage, args.out), indent=2))
    elif args.command == "reversible-candidate-lineage-carbon":
        from .candidate_lineage import attach_candidate_carbon_mapping
        print(json.dumps(attach_candidate_carbon_mapping(args.lineage, args.mapping, args.out), indent=2))
    elif args.command == "map-identity-set-upstream":
        print(json.dumps(map_identity_set_upstream(args.source, args.out), indent=2))
    elif args.command == "identity-set-core-bridges":
        print(json.dumps(build_identity_set_core_bridges(args.source, args.network, args.out), indent=2))
    elif args.command == "extract-cannabisdb-table":
        print(json.dumps(extract_terpedia_table(args.xml, args.out, args.report), indent=2))


if __name__ == "__main__":
    main()
