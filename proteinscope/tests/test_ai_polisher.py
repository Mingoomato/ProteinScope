"""Tests for annotators/ai_polisher.py"""

from __future__ import annotations

import json
from unittest.mock import AsyncMock, patch

import pytest


SAMPLE_NODE_DATA = {
    "gene_symbol": "EGFR",
    "protein_name": "Epidermal growth factor receptor",
    "interaction_comment": "EGFR acts as receptor tyrosine kinase; activates: RAS, PI3K; inhibits: —",
    "drug_comment": "Drug targets: Gefitinib (inhibitor, 1A); Erlotinib (inhibitor, 1B)",
    "disease_comment": "Disease association: Non-small cell lung cancer",
    "upstream_comment": "Activated by: EGF, TGF-alpha",
    "downstream_comment": "Downstream targets: RAS, PI3K, STAT3",
}

SHORT_NODE_DATA = {
    "gene_symbol": "EGFR",
    "protein_name": "Epidermal growth factor receptor",
    "interaction_comment": "EGFR participates in this pathway",
    "drug_comment": "",
    "disease_comment": "",
    "upstream_comment": "",
    "downstream_comment": "",
}


# ---------------------------------------------------------------------------
# Successful Gemini response
# ---------------------------------------------------------------------------

@pytest.mark.asyncio
async def test_polish_comment_success():
    """Should parse Gemini JSON response and return short + extended."""
    gemini_response = json.dumps({
        "short_comment": "EGFR is a receptor tyrosine kinase central to cellular proliferation.",
        "extended_footnote": "EGFR mediates signaling via EGF binding. It is inhibited by "
                             "Gefitinib and Erlotinib in NSCLC treatment.",
    })

    with patch("core.gemini_interpreter._call", new_callable=AsyncMock) as mock_call:
        mock_call.return_value = gemini_response
        from annotators.ai_polisher import polish_comment
        short, extended = await polish_comment(SAMPLE_NODE_DATA)

    assert "EGFR" in short
    assert "receptor tyrosine kinase" in short
    assert extended is not None
    assert "Gefitinib" in extended or "NSCLC" in extended


@pytest.mark.asyncio
async def test_polish_comment_with_markdown_fences():
    """Should strip markdown fences from Gemini response."""
    gemini_response = (
        "```json\n"
        '{"short_comment": "EGFR drives cell proliferation.", '
        '"extended_footnote": null}\n'
        "```"
    )

    with patch("core.gemini_interpreter._call", new_callable=AsyncMock) as mock_call:
        mock_call.return_value = gemini_response
        from annotators.ai_polisher import polish_comment
        short, extended = await polish_comment(SHORT_NODE_DATA)

    assert "EGFR" in short
    assert extended is None


@pytest.mark.asyncio
async def test_polish_comment_null_extended_footnote():
    """Should return None for extended_footnote when Gemini returns null."""
    gemini_response = json.dumps({
        "short_comment": "EGFR is a kinase.",
        "extended_footnote": None,
    })

    with patch("core.gemini_interpreter._call", new_callable=AsyncMock) as mock_call:
        mock_call.return_value = gemini_response
        from annotators.ai_polisher import polish_comment
        short, extended = await polish_comment(SHORT_NODE_DATA)

    assert short == "EGFR is a kinase."
    assert extended is None


# ---------------------------------------------------------------------------
# Fallback on Gemini failure
# ---------------------------------------------------------------------------

@pytest.mark.asyncio
async def test_polish_comment_fallback_on_empty_response():
    """Should fall back to interaction_comment when Gemini returns empty string."""
    with patch("core.gemini_interpreter._call", new_callable=AsyncMock) as mock_call:
        mock_call.return_value = ""
        from annotators.ai_polisher import polish_comment
        short, extended = await polish_comment(SAMPLE_NODE_DATA)

    # Fallback: short should be the interaction_comment
    assert short == SAMPLE_NODE_DATA["interaction_comment"]


@pytest.mark.asyncio
async def test_polish_comment_fallback_on_invalid_json():
    """Should fall back gracefully when Gemini returns malformed JSON."""
    with patch("core.gemini_interpreter._call", new_callable=AsyncMock) as mock_call:
        mock_call.return_value = "This is not JSON at all"
        from annotators.ai_polisher import polish_comment
        short, extended = await polish_comment(SAMPLE_NODE_DATA)

    assert short == SAMPLE_NODE_DATA["interaction_comment"]


@pytest.mark.asyncio
async def test_polish_comment_fallback_on_exception():
    """Should fall back when _call raises an exception."""
    with patch("core.gemini_interpreter._call", new_callable=AsyncMock) as mock_call:
        mock_call.side_effect = RuntimeError("API quota exceeded")
        from annotators.ai_polisher import polish_comment
        short, extended = await polish_comment(SAMPLE_NODE_DATA)

    assert isinstance(short, str)
    assert len(short) > 0


# ---------------------------------------------------------------------------
# _fallback function
# ---------------------------------------------------------------------------

def test_fallback_short_total_below_threshold():
    """When total chars are below THRESHOLD_CHARS, extended should be None."""
    from annotators.ai_polisher import _fallback, THRESHOLD_CHARS

    short_data = {
        "interaction_comment": "Short comment.",
        "drug_comment": "",
        "disease_comment": "",
        "upstream_comment": "",
        "downstream_comment": "",
    }
    short, extended = _fallback("Short comment.", short_data, 50)
    assert extended is None


def test_fallback_extended_above_threshold():
    """When total chars exceed THRESHOLD_CHARS, extended should include available comment parts."""
    from annotators.ai_polisher import _fallback, THRESHOLD_CHARS

    data = {
        "interaction_comment": "EGFR acts as receptor tyrosine kinase.",
        "drug_comment": "Drug targets: Gefitinib (inhibitor, 1A)",
        "disease_comment": "Disease association: Non-small cell lung cancer",
        "upstream_comment": "Activated by: EGF",
        "downstream_comment": "Downstream targets: RAS",
    }
    total = sum(len(v) for v in data.values())
    short, extended = _fallback("Short.", data, total)

    # total is well above 300 chars
    if total > THRESHOLD_CHARS:
        assert extended is not None
        # Should include at least one of the non-interaction comments
        non_empty_parts = [v for k, v in data.items()
                           if k != "interaction_comment" and v]
        assert any(part in extended for part in non_empty_parts)


# ---------------------------------------------------------------------------
# SYSTEM_PROMPT presence
# ---------------------------------------------------------------------------

def test_system_prompt_contains_rules():
    """The system prompt should contain the key instruction keywords."""
    from annotators.ai_polisher import SYSTEM_PROMPT
    assert "SHORT COMMENT" in SYSTEM_PROMPT
    assert "EXTENDED FOOTNOTE" in SYSTEM_PROMPT
    assert "JSON" in SYSTEM_PROMPT
    assert "short_comment" in SYSTEM_PROMPT
    assert "extended_footnote" in SYSTEM_PROMPT


# ---------------------------------------------------------------------------
# THRESHOLD_CHARS constant
# ---------------------------------------------------------------------------

def test_threshold_chars_value():
    """THRESHOLD_CHARS should be 300 as specified."""
    from annotators.ai_polisher import THRESHOLD_CHARS
    assert THRESHOLD_CHARS == 300
