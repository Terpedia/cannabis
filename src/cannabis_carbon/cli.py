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
    p_complete.add_argument("--out", type=Path, default=Path("data/reports/completeness.json"))
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
        result = compute_completeness(args.network, args.compounds, args.mapping)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
