"""Tests for sequence_aligner.py"""

from __future__ import annotations

import pytest
from analyzers.sequence_aligner import (
    align_pairwise,
    align_at_annotated_sites,
    build_msa,
    AlignmentResult,
    MSAResult,
)


def test_identical_sequences():
    seq = (
        "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQTLSEQVQEELLSSQVTQELR"
        "ALMDETMKELKAYKSELEEQLTPVAEETRARLSKELQAAQARLGADVLASHGRLVQYRGEVQAMLGQSTEELRVRLASHLR"
        "KLRKRLLRDA"
    )
    result = align_pairwise(seq, seq)
    assert result is not None
    assert result.identity_pct == pytest.approx(100.0, abs=0.1)


def test_different_sequences():
    seq_a = "MVKVLWALVTFLAGA"
    seq_b = "MKVLWALLVTFLAGC"
    result = align_pairwise(seq_a, seq_b)
    assert result is not None
    assert 0 < result.identity_pct < 100
    assert result.gap_pct >= 0


def test_empty_sequence_returns_none():
    result = align_pairwise("", "MKVLW")
    assert result is None


def test_both_empty_returns_none():
    result = align_pairwise("", "")
    assert result is None


def test_alignment_result_fields():
    seq_a = "MKVLWAALLV"
    seq_b = "MKVLWAALLV"
    result = align_pairwise(seq_a, seq_b)
    assert result is not None
    assert isinstance(result, AlignmentResult)
    assert 0.0 <= result.identity_pct <= 100.0
    assert 0.0 <= result.similarity_pct <= 100.0
    assert 0.0 <= result.gap_pct <= 100.0
    assert isinstance(result.aligned_query, str)
    assert isinstance(result.aligned_target, str)
    assert isinstance(result.score, float)


def test_similarity_gte_identity():
    """Similarity should be >= identity since it includes conservative substitutions."""
    seq_a = "MKVLWAALLVTFLAG"
    seq_b = "MKVLWALLLTFLAGC"
    result = align_pairwise(seq_a, seq_b)
    assert result is not None
    assert result.similarity_pct >= result.identity_pct


def test_site_conservation():
    from analyzers.sequence_aligner import align_at_annotated_sites
    seq_a = "MKVLWAALLV"
    seq_b = "MKVLWAALLV"
    sites = [{"start": 0, "end": 3, "description": "binding site", "sequence_fragment": "MKV"}]
    result = align_at_annotated_sites(seq_a, seq_b, sites)
    assert len(result) == 3  # positions 0, 1, 2
    for entry in result:
        assert entry["conserved"] is True


def test_site_conservation_identical_single_residue():
    from analyzers.sequence_aligner import align_at_annotated_sites
    seq_a = "MKVLWAALLV"
    seq_b = "MKVLWAALLV"
    sites = [{"start": 0, "end": 1, "description": "start_met"}]
    result = align_at_annotated_sites(seq_a, seq_b, sites)
    assert len(result) == 1
    assert result[0]["conserved"] is True
    assert result[0]["residue_a"] == "M"
    assert result[0]["residue_b"] == "M"


def test_site_conservation_empty_sites():
    from analyzers.sequence_aligner import align_at_annotated_sites
    result = align_at_annotated_sites("MKVLW", "MKVLW", [])
    assert result == []


def test_msa_returns_result():
    seqs = {
        "seq_a": "MKVLWAALLVTFLAG",
        "seq_b": "MKVLWAALLVTFLAG",
        "seq_c": "MKVLWAALLATFLAG",
    }
    result = build_msa(seqs)
    assert result is not None
    assert result.n_sequences == 3
    assert len(result.conservation_scores) > 0


def test_msa_conservation_scores_range():
    """All conservation scores should be in [0.0, 1.0]."""
    seqs = {
        "alpha": "MKVLWAALLVTFLAG",
        "beta":  "MKVLWAALLVTFLAG",
        "gamma": "MKVLWAALLATFLAG",
    }
    result = build_msa(seqs)
    assert result is not None
    for score in result.conservation_scores:
        assert 0.0 <= score <= 1.0


def test_msa_identical_sequences_fully_conserved():
    """Identical sequences should produce conservation scores close to 1.0."""
    seqs = {
        "x": "MKVLWAALLV",
        "y": "MKVLWAALLV",
    }
    result = build_msa(seqs)
    assert result is not None
    for score in result.conservation_scores:
        assert score == pytest.approx(1.0, abs=0.01)


def test_msa_single_sequence_returns_none():
    result = build_msa({"only": "MKVLW"})
    assert result is None


def test_msa_empty_returns_none():
    result = build_msa({})
    assert result is None


def test_msa_result_type():
    seqs = {"a": "MKVLWAALLV", "b": "MKVLWAALLV"}
    result = build_msa(seqs)
    assert isinstance(result, MSAResult)
    assert isinstance(result.aligned_sequences, dict)
    assert isinstance(result.conservation_scores, list)
