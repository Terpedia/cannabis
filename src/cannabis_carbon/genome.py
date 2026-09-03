"""Attach reference-proteome sequence evidence to Terpedia candidate hypotheses."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


def _fasta(path: Path) -> dict[str, str]:
    sequences, accession, chunks = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if accession:
                sequences[accession] = "".join(chunks)
            match = re.search(r"(?:sp|tr)\|([A-Z0-9]+)\|", line)
            accession = match.group(1) if match else line[1:].split()[0].split("|")[-1]
            chunks = []
        else:
            chunks.append(line.strip())
    if accession:
        sequences[accession] = "".join(chunks)
    return sequences


def build_genome_search(queue_path: Path, fasta_path: Path, output: Path) -> dict:
    queue = json.loads(queue_path.read_text())
    sequences = _fasta(fasta_path)
    candidates = {}
    for item in queue["items"]:
        for protein in item.get("candidate_proteins", []):
            protein_id = protein.get("proteinId")
            accession = protein.get("accession")
            if not protein_id:
                continue
            sequence = sequences.get(accession, "")
            candidates[protein_id] = {**protein, "sequence_search": {"method": "reference-proteome-membership", "proteome": "UP000583929", "fasta_source": str(fasta_path), "sequence_present": bool(sequence), "length": len(sequence) or None, "sha256": hashlib.sha256(sequence.encode()).hexdigest() if sequence else None, "claim": "Sequence presence and existing annotation support a candidate protein; they do not establish reaction specificity or in-vivo activity."}}
    report = {"schema": "cannabis-carbon.genome-candidate-search.v1", "source_queue": str(queue_path), "proteome": {"accession": "UP000583929", "fasta": str(fasta_path), "protein_count": len(sequences), "fasta_sha256": hashlib.sha256(fasta_path.read_bytes()).hexdigest()}, "method": "reference-proteome-membership", "candidate_protein_count": len(candidates), "candidate_proteins_with_sequence": sum(bool(p["sequence_search"]["sequence_present"]) for p in candidates.values()), "candidate_proteins_without_sequence": sum(not p["sequence_search"]["sequence_present"] for p in candidates.values()), "candidate_proteins": list(candidates.values()), "hypotheses": queue["items"], "claim_boundary": "This report records genome/proteome membership and annotation evidence. It is not a BLAST, Diamond homology, motif, domain, or biochemical specificity result; candidates must not be promoted automatically."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {"proteome_protein_count": len(sequences), "candidate_protein_count": len(candidates), "candidate_proteins_with_sequence": report["candidate_proteins_with_sequence"], "candidate_proteins_without_sequence": report["candidate_proteins_without_sequence"]}
