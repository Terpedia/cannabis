"""Rank genome-derived enzyme candidates without over-claiming function.

The expensive searches (DIAMOND/HMMER, motif and domain annotation) are kept
outside this small library. Their outputs are normalized into a JSON record and
combined here so ranking is deterministic and auditable.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict


@dataclass(frozen=True)
class CandidateEvidence:
    protein_id: str
    reaction_id: str
    identity: float | None = None
    coverage: float | None = None
    profile_score: float | None = None
    catalytic_motif: bool = False
    complete_domains: bool = False
    localization_support: bool = False
    expression_support: bool = False
    source_urls: tuple[str, ...] = ()


def _bounded(value: float | None) -> float:
    return max(0.0, min(1.0, value or 0.0))


def rank_candidate(candidate: CandidateEvidence) -> dict:
    """Return a weighted score and evidence-qualified status.

    Identity and coverage are intentionally capped: high homology cannot
    establish substrate specificity. Only direct biochemical evidence can be
    labeled ``biochemically_supported`` by downstream curation.
    """
    score = (
        0.30 * _bounded(candidate.identity)
        + 0.20 * _bounded(candidate.coverage)
        + 0.25 * _bounded(candidate.profile_score)
        + 0.10 * candidate.catalytic_motif
        + 0.07 * candidate.complete_domains
        + 0.04 * candidate.localization_support
        + 0.04 * candidate.expression_support
    )
    status = "homology_candidate" if score >= 0.45 else "weak_candidate"
    return {**asdict(candidate), "score": round(score, 6), "status": status}
