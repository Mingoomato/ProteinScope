"""Tests for fetchers/kegg.py using pytest-httpx."""

from __future__ import annotations

import pytest

try:
    from pytest_httpx import HTTPXMock
    HAS_HTTPX_MOCK = True
except ImportError:
    HAS_HTTPX_MOCK = False

from fetchers.kegg import (
    get_pathways_for_gene,
    get_disease_pathways_for_gene,
    get_pathway_detail,
    get_reactions_for_pathway,
    get_reaction_detail,
    parse_kegg_flat,
    fetch_all_kegg,
)

pytestmark = pytest.mark.skipif(
    not HAS_HTTPX_MOCK,
    reason="pytest-httpx not installed",
)

KEGG_BASE = "https://rest.kegg.jp"


@pytest.mark.asyncio
async def test_get_pathways_for_gene_ldha(httpx_mock: "HTTPXMock"):
    """NCBI gene 3939 (LDHA) should map to glycolysis pathway."""
    httpx_mock.add_response(
        url=f"{KEGG_BASE}/link/pathway/hsa:3939",
        text="hsa:3939\tpath:hsa00010\nhsa:3939\tpath:hsa01200\n",
    )
    result = await get_pathways_for_gene("3939")
    assert "path:hsa00010" in result
    assert len(result) == 2


@pytest.mark.asyncio
async def test_get_disease_pathways_for_gene(httpx_mock: "HTTPXMock"):
    """CFTR gene should be linked to cystic fibrosis disease entry."""
    httpx_mock.add_response(
        url=f"{KEGG_BASE}/link/disease/hsa:1080",
        text="hsa:1080\tds:H00218\n",
    )
    result = await get_disease_pathways_for_gene("1080")
    assert "ds:H00218" in result


@pytest.mark.asyncio
async def test_parse_kegg_flat():
    """parse_kegg_flat should correctly parse KEGG flat-file format."""
    flat_text = (
        "ENTRY       hsa00010          Pathway\n"
        "NAME        Glycolysis / Gluconeogenesis\n"
        "ORGANISM    Homo sapiens (human) [GN:hsa]\n"
        "///"
    )
    result = parse_kegg_flat(flat_text)
    assert "NAME" in result
    assert "Glycolysis" in result["NAME"]


@pytest.mark.asyncio
async def test_get_pathway_detail(httpx_mock: "HTTPXMock"):
    """Should return flat-file text for a pathway."""
    httpx_mock.add_response(
        url=f"{KEGG_BASE}/get/hsa00010",
        text="ENTRY       hsa00010\nNAME        Glycolysis\n///",
    )
    result = await get_pathway_detail("path:hsa00010")
    assert "Glycolysis" in result


@pytest.mark.asyncio
async def test_get_reaction_detail(httpx_mock: "HTTPXMock"):
    """Should return parsed reaction entry."""
    httpx_mock.add_response(
        url=f"{KEGG_BASE}/get/R00299",
        text=(
            "ENTRY       R00299\n"
            "NAME        ATP:pyruvate 2-O-phosphotransferase\n"
            "EQUATION    ATP + Pyruvate <=> ADP + PEP\n"
            "///"
        ),
    )
    result = await get_reaction_detail("rn:R00299")
    assert "EQUATION" in result
    assert "ATP" in result["EQUATION"]


@pytest.mark.asyncio
async def test_fetch_all_kegg_empty_gene_id():
    """Empty NCBI gene ID should return empty structure without API calls."""
    result = await fetch_all_kegg("")
    assert result["pathway_ids"] == []
    assert result["disease_pathway_ids"] == []


@pytest.mark.asyncio
async def test_fetch_all_kegg_integration(httpx_mock: "HTTPXMock"):
    """Full fetch_all_kegg for LDHA (ncbi_gene_id=3939)."""
    httpx_mock.add_response(
        url=f"{KEGG_BASE}/link/pathway/hsa:3939",
        text="hsa:3939\tpath:hsa00010\n",
    )
    httpx_mock.add_response(
        url=f"{KEGG_BASE}/link/disease/hsa:3939",
        text="",
    )
    httpx_mock.add_response(
        url=f"{KEGG_BASE}/get/hsa00010",
        text="ENTRY       hsa00010\nNAME        Glycolysis\n///",
    )
    result = await fetch_all_kegg("3939")
    assert "path:hsa00010" in result["pathway_ids"]
    assert "hsa00010" in result["pathway_details"]
