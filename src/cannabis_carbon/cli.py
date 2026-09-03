from __future__ import annotations

import argparse
import hashlib
import json
import ssl
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem

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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p_download = sub.add_parser("download")
    p_download.add_argument("--out", type=Path, default=Path("data/raw"))
    p_download.add_argument("--insecure-download", action="store_true", help="Allow the current CannabisDB expired TLS certificate")
    p_inspect = sub.add_parser("inspect-sdf")
    p_inspect.add_argument("path", type=Path)
    args = parser.parse_args()
    if args.command == "download":
        download(args.out, args.insecure_download)
    elif args.command == "inspect-sdf":
        inspect_sdf(args.path)
