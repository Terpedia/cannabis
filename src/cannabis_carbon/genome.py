"""Attach reference-proteome sequence evidence to Terpedia candidate hypotheses."""

from __future__ import annotations

import hashlib
import json
import re
import csv
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


def build_genome_search(queue_path: Path, fasta_path: Path, output: Path, diamond_hits_path: Path | None = None, reference_tsv_path: Path | None = None) -> dict:
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
    reference_ec = {}
    if reference_tsv_path and reference_tsv_path.exists():
        with reference_tsv_path.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                reference_ec[row["Entry"]] = [ec for ec in row.get("EC number", "").split(";") if ec]
    if diamond_hits_path and diamond_hits_path.exists():
        hits_by_accession = {}
        for line in diamond_hits_path.read_text().splitlines():
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            query = fields[0].split("|")[1] if "|" in fields[0] else fields[0]
            subject_parts = fields[1].split("|")
            subject = subject_parts[1] if len(subject_parts) > 2 else fields[1]
            hit = {"reference_accession": subject, "reference_ec_numbers": reference_ec.get(subject, []), "identity_percent": float(fields[2]), "alignment_length": int(fields[3]), "query_length": int(fields[4]), "reference_length": int(fields[5]), "evalue": float(fields[6]), "bitscore": float(fields[7]), "query_coverage_percent": float(fields[8]), "reference_coverage_percent": float(fields[9]), "source": str(diamond_hits_path)}
            hits_by_accession.setdefault(query, []).append(hit)
        for protein in candidates.values():
            accession = protein.get("accession")
            hits = sorted(hits_by_accession.get(accession, []), key=lambda h: (h["evalue"], -h["bitscore"]))[:10]
            if hits:
                protein["diamond_search"] = {"method": "DIAMOND blastp", "hits": hits, "hit_count_reported": len(hits), "thresholds": {"evalue_max": 1e-5, "identity_min_for_strong_candidate": 30.0, "query_coverage_min_for_strong_candidate": 50.0}, "strong_candidate_hit": any(h["identity_percent"] >= 30.0 and h["query_coverage_percent"] >= 50.0 and h["evalue"] <= 1e-5 for h in hits), "claim": "Sequence similarity is candidate evidence only; it does not establish substrate specificity or in-vivo activity."}
            elif protein.get("specialized_search"):
                protein["diamond_search"] = protein["specialized_search"]
    for protein in candidates.values():
        if not protein.get("diamond_search") and protein.get("specialized_search"):
            protein["diamond_search"] = protein["specialized_search"]
    diamond_present = [p for p in candidates.values() if p.get("diamond_search")]
    report = {"schema": "cannabis-carbon.genome-candidate-search.v2", "source_queue": str(queue_path), "proteome": {"accession": "UP000583929", "fasta": str(fasta_path), "protein_count": len(sequences), "fasta_sha256": hashlib.sha256(fasta_path.read_bytes()).hexdigest()}, "method": "reference-proteome-membership plus DIAMOND blastp when hit file is supplied", "reference_database": {"tsv": str(reference_tsv_path) if reference_tsv_path else None, "reference_sequence_count": sum(1 for _ in reference_ec)}, "diamond_hits": str(diamond_hits_path) if diamond_hits_path else None, "candidate_protein_count": len(candidates), "candidate_proteins_with_sequence": sum(bool(p["sequence_search"]["sequence_present"]) for p in candidates.values()), "candidate_proteins_without_sequence": sum(not p["sequence_search"]["sequence_present"] for p in candidates.values()), "candidate_proteins_with_diamond_hits": sum(bool(p.get("diamond_search", {}).get("hits")) for p in candidates.values()), "candidate_proteins_meeting_strong_hit_threshold": sum(bool(p.get("diamond_search", {}).get("strong_candidate_hit")) for p in candidates.values()), "candidate_proteins": list(candidates.values()), "hypotheses": queue["items"], "claim_boundary": "This report records proteome membership, annotation, and optional sequence-similarity evidence. Candidates must not be promoted automatically to confirmed enzymes."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    return {"proteome_protein_count": len(sequences), "candidate_protein_count": len(candidates), "candidate_proteins_with_sequence": report["candidate_proteins_with_sequence"], "candidate_proteins_without_sequence": report["candidate_proteins_without_sequence"]}
