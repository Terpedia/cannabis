from cannabis_carbon.candidates import CandidateEvidence, rank_candidate


def test_candidate_is_not_confirmed_by_homology():
    result = rank_candidate(CandidateEvidence(
        protein_id="p1", reaction_id="r1", identity=0.9, coverage=0.95,
        profile_score=0.9, catalytic_motif=True, complete_domains=True,
    ))
    assert result["status"] == "homology_candidate"
    assert result["status"] != "biochemically_supported"


def test_weak_candidate_is_retained():
    result = rank_candidate(CandidateEvidence(protein_id="p2", reaction_id="r2"))
    assert result["status"] == "weak_candidate"
