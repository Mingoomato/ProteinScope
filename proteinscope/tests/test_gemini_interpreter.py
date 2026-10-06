"""Tests for core/gemini_interpreter.py — all Gemini calls are monkeypatched."""

from __future__ import annotations

import pytest
from unittest.mock import AsyncMock, MagicMock, patch


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# resolve_nl_query
# ---------------------------------------------------------------------------

class TestResolveNlQuery:
    @pytest.mark.asyncio
    async def test_warfarin_enzyme_resolves_to_cyp2c9(self):
        with patch("core.gemini_interpreter._call_fast", new_callable=AsyncMock, return_value="CYP2C9"):
            from core.gemini_interpreter import resolve_nl_query
            result = await resolve_nl_query("the enzyme that breaks down warfarin")
        assert result == "CYP2C9"

    @pytest.mark.asyncio
    async def test_cftr_description_resolves(self):
        with patch("core.gemini_interpreter._call_fast", new_callable=AsyncMock, return_value="CFTR"):
            from core.gemini_interpreter import resolve_nl_query
            result = await resolve_nl_query("cystic fibrosis channel protein")
        assert result == "CFTR"

    @pytest.mark.asyncio
    async def test_unknown_returns_none(self):
        with patch("core.gemini_interpreter._call_fast", new_callable=AsyncMock, return_value="UNKNOWN"):
            from core.gemini_interpreter import resolve_nl_query
            result = await resolve_nl_query("something completely unrecognizable xyz123")
        assert result is None

    @pytest.mark.asyncio
    async def test_empty_response_returns_none(self):
        with patch("core.gemini_interpreter._call_fast", new_callable=AsyncMock, return_value=""):
            from core.gemini_interpreter import resolve_nl_query
            result = await resolve_nl_query("any query")
        assert result is None

    @pytest.mark.asyncio
    async def test_no_api_key_returns_none(self):
        with patch("core.gemini_interpreter._call_fast", new_callable=AsyncMock, return_value=""):
            from core.gemini_interpreter import resolve_nl_query
            result = await resolve_nl_query("EGFR receptor")
        assert result is None

    @pytest.mark.asyncio
    async def test_strips_punctuation_from_response(self):
        with patch("core.gemini_interpreter._call_fast", new_callable=AsyncMock, return_value="EGFR."):
            from core.gemini_interpreter import resolve_nl_query
            result = await resolve_nl_query("HER1 cancer receptor")
        assert result == "EGFR"


# ---------------------------------------------------------------------------
# generate_report_summary
# ---------------------------------------------------------------------------

class TestGenerateReportSummary:
    def _make_record(self):
        """Build a minimal ProteinRecord-like object for testing."""
        record = MagicMock()
        record.protein_name = "Cytochrome P450 2C9"
        record.gene_name = "CYP2C9"
        record.organism = "Homo sapiens"
        record.function_description = "Oxidizes small foreign organic molecules."
        record.clinical_variants = []
        record.drug_interactions = []
        record.metabolic_pathways = []
        record.signaling_pathways = []
        record.pgx_variants = []
        return record

    @pytest.mark.asyncio
    async def test_returns_summary_string(self):
        expected = "CYP2C9 is a key drug-metabolizing enzyme..."
        with patch("core.gemini_interpreter._call", new_callable=AsyncMock, return_value=expected):
            from core.gemini_interpreter import generate_report_summary
            result = await generate_report_summary(self._make_record())
        assert result == expected

    @pytest.mark.asyncio
    async def test_no_key_returns_empty_string(self):
        with patch("core.gemini_interpreter._call", new_callable=AsyncMock, return_value=""):
            from core.gemini_interpreter import generate_report_summary
            result = await generate_report_summary(self._make_record())
        assert result == ""

    @pytest.mark.asyncio
    async def test_api_failure_returns_empty_string(self):
        with patch("core.gemini_interpreter._call", new_callable=AsyncMock, return_value=""):
            from core.gemini_interpreter import generate_report_summary
            result = await generate_report_summary(self._make_record())
        assert result == ""


# ---------------------------------------------------------------------------
# generate_section_insight
# ---------------------------------------------------------------------------

class TestGenerateSectionInsight:
    @pytest.mark.asyncio
    async def test_returns_insight_string(self):
        expected = "CYP2C9 metabolizes ~15% of all clinical drugs."
        with patch("core.gemini_interpreter._call", new_callable=AsyncMock, return_value=expected):
            from core.gemini_interpreter import generate_section_insight
            result = await generate_section_insight(
                "Drug Interactions",
                "CYP2C9 interacts with warfarin, celecoxib, and phenytoin."
            )
        assert result == expected

    @pytest.mark.asyncio
    async def test_no_key_returns_empty_string(self):
        with patch("core.gemini_interpreter._call", new_callable=AsyncMock, return_value=""):
            from core.gemini_interpreter import generate_section_insight
            result = await generate_section_insight("Any Section", "any data")
        assert result == ""
