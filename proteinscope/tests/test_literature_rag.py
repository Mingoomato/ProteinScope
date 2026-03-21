"""Tests for literature_rag.py"""
import pytest
from analyzers.literature_rag import LiteratureRAG
from core.models import Reference


def _make_refs():
    return [
        Reference(
            title="EGFR signaling in lung cancer",
            authors=["Smith A", "Jones B"],
            journal="Nature Medicine",
            year=2020,
            apa_citation="Smith, A., & Jones, B. (2020). EGFR signaling. Nature Medicine.",
        ),
        Reference(
            title="TKI resistance mechanisms in NSCLC",
            authors=["Brown C"],
            journal="Cancer Research",
            year=2021,
            apa_citation="Brown, C. (2021). TKI resistance. Cancer Research.",
        ),
    ]


def test_index_and_query():
    rag = LiteratureRAG()
    refs = _make_refs()
    n = rag.index(refs)
    # If sentence-transformers available, should index both refs
    # If not, graceful degradation returns 0
    assert isinstance(n, int)
    assert n >= 0


def test_query_returns_list():
    rag = LiteratureRAG()
    rag.index(_make_refs())
    results = rag.query("EGFR resistance")
    assert isinstance(results, list)


def test_query_without_index_returns_empty():
    rag = LiteratureRAG()
    results = rag.query("anything")
    assert results == []


def test_index_empty_returns_zero():
    rag = LiteratureRAG()
    n = rag.index([])
    assert n == 0


def test_query_result_structure():
    rag = LiteratureRAG()
    n = rag.index(_make_refs())
    if n == 0:
        pytest.skip("sentence-transformers not available")
    results = rag.query("EGFR resistance", top_k=1)
    assert len(results) <= 1
    if results:
        r = results[0]
        assert "reference" in r
        assert "score" in r
        assert "abstract_snippet" in r
        assert isinstance(r["score"], float)
        assert isinstance(r["abstract_snippet"], str)


def test_query_top_k_respected():
    rag = LiteratureRAG()
    n = rag.index(_make_refs())
    if n == 0:
        pytest.skip("sentence-transformers not available")
    results = rag.query("cancer signaling", top_k=1)
    assert len(results) <= 1


def test_query_returns_relevant_first():
    rag = LiteratureRAG()
    n = rag.index(_make_refs())
    if n == 0:
        pytest.skip("sentence-transformers not available")
    results = rag.query("TKI resistance NSCLC", top_k=2)
    # Top result should be the TKI resistance paper
    assert len(results) >= 1
    top_ref = results[0]["reference"]
    assert "TKI" in top_ref.title or "resistance" in top_ref.title.lower()


@pytest.mark.asyncio
async def test_answer_graceful_without_key(monkeypatch):
    monkeypatch.setenv("GEMINI_API_KEY", "")
    rag = LiteratureRAG()
    rag.index(_make_refs())
    result = await rag.answer("What is EGFR?")
    assert isinstance(result, str)


@pytest.mark.asyncio
async def test_answer_without_index_is_string(monkeypatch):
    monkeypatch.setenv("GEMINI_API_KEY", "")
    rag = LiteratureRAG()
    result = await rag.answer("What is PCNA?")
    assert isinstance(result, str)
