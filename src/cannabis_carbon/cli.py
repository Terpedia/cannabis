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
from .hypotheses import build_candidate_queue
from .crosswalk import build_crosswalk
from .balance import audit_balances
from .lineage import build_carbon_lineage
from .networkdb import build_networkdb
from .genome import build_genome_search
from .inventory import build_specialty_inventory
from .test_hypotheses import build_test_hypotheses

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
    p_complete.add_argument("--out", type=Path, default=Path("data/reports/completeness.json"))
    p_queue = sub.add_parser("candidate-queue")
    p_queue.add_argument("source", type=Path)
    p_queue.add_argument("--out", type=Path, default=Path("data/reports/candidate-work-queue.json"))
    p_crosswalk = sub.add_parser("crosswalk")
    p_crosswalk.add_argument("cannabisdb_sdf", type=Path)
    p_crosswalk.add_argument("terpedia_network", type=Path)
    p_crosswalk.add_argument("--out", type=Path, default=Path("data/reports/identity-crosswalk.json"))
    p_balance = sub.add_parser("balance-audit")
    p_balance.add_argument("network", type=Path)
    p_balance.add_argument("--out", type=Path, default=Path("data/reports/phase1-balance-audit.json"))
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
    p_networkdb.add_argument("--out", type=Path, default=Path("docs/data/networkdb.json"))
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
    p_inventory.add_argument("--out", type=Path, default=Path("data/reports/named-specialty-inventory.json"))
    p_tests = sub.add_parser("test-hypotheses")
    p_tests.add_argument("queue", type=Path)
    p_tests.add_argument("network", type=Path)
    p_tests.add_argument("--lineage", type=Path, default=Path("data/reports/carbon-lineage.json"))
    p_tests.add_argument("--compounds", type=Path, default=Path("docs/data/compounds.json"))
    p_tests.add_argument("--out", type=Path, default=Path("data/reports/testable-hypotheses.json"))
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
        graph = cytoscape_elements(load_network(args.source))
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(graph, separators=(",", ":")) + "\n")
        print(json.dumps(graph["stats"], indent=2))
    elif args.command == "map-reactions":
        print(json.dumps(build_reaction_report(args.source, args.out), indent=2))
    elif args.command == "completeness":
        result = compute_completeness(args.network, args.compounds, args.mapping, args.crosswalk, args.lineage, args.atom_audit)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
    elif args.command == "candidate-queue":
        print(json.dumps(build_candidate_queue(args.source, args.out), indent=2))
    elif args.command == "crosswalk":
        print(json.dumps(build_crosswalk(args.cannabisdb_sdf, args.terpedia_network, args.out), indent=2))
    elif args.command == "balance-audit":
        print(json.dumps(audit_balances(args.network, args.out), indent=2))
    elif args.command == "carbon-lineage":
        print(json.dumps(build_carbon_lineage(args.network, args.mapping, args.crosswalk, args.compounds, args.out, args.directions), indent=2))
    elif args.command == "carbon-atom-audit":
        from .lineage import build_carbon_atom_audit
        print(json.dumps(build_carbon_atom_audit(args.network, args.lineage, args.crosswalk, args.compounds, args.out), indent=2))
    elif args.command == "networkdb":
        print(json.dumps(build_networkdb(args.network, args.compounds, args.crosswalk, args.out, args.hypotheses, args.genome_search, args.genome_fasta, args.mapping, args.lineage, args.atom_audit), indent=2))
    elif args.command == "genome-search":
        print(json.dumps(build_genome_search(args.queue, args.fasta, args.out, args.diamond_hits, args.reference_tsv), indent=2))
    elif args.command == "specialty-inventory":
        print(json.dumps(build_specialty_inventory(args.compounds, args.crosswalk, args.network, args.out), indent=2))
    elif args.command == "test-hypotheses":
        print(json.dumps(build_test_hypotheses(args.queue, args.network, args.out, args.lineage, args.compounds), indent=2))


if __name__ == "__main__":
    main()
