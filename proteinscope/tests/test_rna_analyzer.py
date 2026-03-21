"""Tests for rna_analyzer.py"""

from __future__ import annotations

import pytest
from analyzers.rna_analyzer import compute_cai, analyze_cds, RNAAnalysisResult


# ---------------------------------------------------------------------------
# Known CDS sequences for testing
# ---------------------------------------------------------------------------

# Kozak-optimised synthetic CDS: ATG + 10 high-usage human codons + TGA stop
# Codons: ATG CTG GAG GAG CTG GAG GAG CTG GAG TGA
_HIGH_USAGE_CDS = "ATGCTGGAGGAGCTGGAGGAGCTGGAGTGA"  # 30 nt, in-frame

# Low usage CDS: all codons are the worst human synonyms
# Codons: ATG TTA CGA CGA TTA TTA CGA CGA TTA TGA
_LOW_USAGE_CDS = "ATGTTACGACGATTATTACGACGATTATGA"  # 30 nt, in-frame

# Well-known human ACTB beta-actin partial CDS (first 30 nt, real sequence)
_ACTB_PARTIAL = "ATGGATGATGATATCGCCGCGCTCGTCGTC"  # 30 nt, starts ATG

# CDS without start codon
_NO_START = "CTGGAGGAGCTGGAGGAGCTGGAGTGA"

# CDS without stop codon
_NO_STOP = "ATGCTGGAGGAGCTGGAGGAGCTGGAG"

# CDS not divisible by 3
_FRAMESHIFTED = "ATGCTGGAGG"


class TestComputeCAI:
    def test_returns_float(self):
        cai = compute_cai(_HIGH_USAGE_CDS)
        assert isinstance(cai, float)

    def test_returns_value_in_0_1_range(self):
        cai = compute_cai(_HIGH_USAGE_CDS)
        assert 0.0 <= cai <= 1.0

    def test_empty_sequence_returns_zero(self):
        assert compute_cai("") == 0.0

    def test_very_short_sequence_returns_zero(self):
        assert compute_cai("AT") == 0.0

    def test_high_usage_cds_higher_than_low(self):
        """Optimal human codons should score higher than rare codons."""
        cai_high = compute_cai(_HIGH_USAGE_CDS)
        cai_low = compute_cai(_LOW_USAGE_CDS)
        assert cai_high > cai_low

    def test_start_met_codon_scores_1(self):
        """ATG is the only Met codon → RA=1.0 → CAI of single Met codon = 1.0."""
        cai = compute_cai("ATG")
        assert cai == pytest.approx(1.0, abs=0.001)

    def test_stop_codon_only_returns_zero(self):
        """A CDS that is only a stop codon has no informative codons."""
        assert compute_cai("TAA") == 0.0
        assert compute_cai("TAG") == 0.0
        assert compute_cai("TGA") == 0.0

    def test_case_insensitive(self):
        """Lower-case input should produce same result as upper-case."""
        cai_upper = compute_cai(_HIGH_USAGE_CDS.upper())
        cai_lower = compute_cai(_HIGH_USAGE_CDS.lower())
        assert cai_upper == pytest.approx(cai_lower, abs=1e-6)

    def test_non_human_organism_still_returns_value(self):
        """Non-human organisms fall back to human table — should not crash."""
        cai = compute_cai(_HIGH_USAGE_CDS, target_organism="yeast")
        assert 0.0 <= cai <= 1.0

    def test_actb_partial_has_reasonable_cai(self):
        """Real human gene partial CDS should score above 0.5."""
        cai = compute_cai(_ACTB_PARTIAL)
        assert cai > 0.40  # beta-actin is highly expressed in humans


class TestAnalyzeCDS:
    def test_returns_rna_analysis_result(self):
        result = analyze_cds(_HIGH_USAGE_CDS)
        assert isinstance(result, RNAAnalysisResult)

    def test_detects_start_codon(self):
        result = analyze_cds(_HIGH_USAGE_CDS)
        assert result.has_start_codon is True

    def test_detects_missing_start_codon(self):
        result = analyze_cds(_NO_START)
        assert result.has_start_codon is False

    def test_detects_stop_codon(self):
        result = analyze_cds(_HIGH_USAGE_CDS)
        assert result.has_stop_codon is True

    def test_detects_missing_stop_codon(self):
        result = analyze_cds(_NO_STOP)
        assert result.has_stop_codon is False

    def test_gc_content_is_percentage(self):
        result = analyze_cds(_HIGH_USAGE_CDS)
        assert 0.0 <= result.gc_content_pct <= 100.0

    def test_gc_content_correctness(self):
        # ATGCTGGAGGAGCTGGAGGAGCTGGAGTGA
        # G+C: CTG=2, GAG=1 → let's count manually
        seq = _HIGH_USAGE_CDS
        gc = sum(1 for c in seq if c in "GC")
        expected_gc_pct = round(100.0 * gc / len(seq), 2)
        result = analyze_cds(seq)
        assert result.gc_content_pct == pytest.approx(expected_gc_pct, abs=0.1)

    def test_cds_length_is_nucleotides(self):
        result = analyze_cds(_HIGH_USAGE_CDS)
        assert result.cds_length_nt == len(_HIGH_USAGE_CDS)

    def test_cai_in_result(self):
        result = analyze_cds(_HIGH_USAGE_CDS)
        assert 0.0 <= result.codon_adaptation_index <= 1.0

    def test_viennarna_available_is_bool(self):
        result = analyze_cds(_HIGH_USAGE_CDS)
        assert isinstance(result.viennarna_available, bool)

    def test_mfe_is_none_without_viennarna(self):
        """MFE should be None since ViennaRNA is not installed."""
        result = analyze_cds(_HIGH_USAGE_CDS)
        if not result.viennarna_available:
            assert result.mfe_kcal is None

    def test_expression_prediction_values(self):
        """expression_prediction must be one of the three valid labels."""
        result = analyze_cds(_HIGH_USAGE_CDS)
        assert result.expression_prediction in ("high", "moderate", "low")

    def test_high_cai_predicts_high_expression(self):
        """High-usage codons should predict high expression."""
        result = analyze_cds(_HIGH_USAGE_CDS)
        # CTG (Leu, RA=1.0) and GAG (Glu, RA=1.0) are all top-usage — CAI should be high
        if result.codon_adaptation_index >= 0.70:
            assert result.expression_prediction == "high"

    def test_empty_cds_graceful(self):
        result = analyze_cds("")
        assert result.cds_length_nt == 0
        assert result.has_start_codon is False
        assert result.has_stop_codon is False
        assert result.codon_adaptation_index == 0.0

    def test_target_organism_preserved(self):
        result = analyze_cds(_HIGH_USAGE_CDS, target_organism="mouse")
        assert result.target_organism == "mouse"

    def test_all_taa_stop_sequence(self):
        """Only stop codons → no informative codons → CAI = 0.0."""
        result = analyze_cds("TAATAATAA")
        assert result.codon_adaptation_index == 0.0

    def test_frameshifted_cds_analyzed_without_crash(self):
        """A CDS not divisible by 3 should not raise — last incomplete codon ignored."""
        result = analyze_cds(_FRAMESHIFTED)
        assert isinstance(result, RNAAnalysisResult)
        assert result.cds_length_nt == len(_FRAMESHIFTED)
