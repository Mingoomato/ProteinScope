"""Tests for fetchers/reactome.py using pytest-httpx."""

from __future__ import annotations

import pytest

try:
    from pytest_httpx import HTTPXMock
    HAS_HTTPX_MOCK = True
except ImportError:
    HAS_HTTPX_MOCK = False

from fetchers.reactome import (
    fetch_pathways,
    fetch_pathway_hierarchy,
    fetch_contained_events,
    fetch_all_reactome,
)

pytestmark = pytest.mark.skipif(
    not HAS_HTTPX_MOCK,
    reason="pytest-httpx not installed",
)

REACTOME_BASE = "https://reactome.org/ContentService"

# Test proteins:
# EGFR  = P00533 (dense signaling pathways)
# CFTR  = P13569 (disease cascade + drug targets)


@pytest.mark.asyncio
async def test_fetch_pathways_egfr(httpx_mock: "HTTPXMock"):
    """P00533 (EGFR) should return signaling pathway entries."""
    httpx_mock.add_response(
        url=f"{REACTOME_BASE}/data/pathways/low/entity/P00533?species=9606",
        json=[
            {"stId": "R-HSA-177929", "displayName": "Signaling by EGFR"},
            {"stId": "R-HSA-1227990", "displayName": "Signaling by FGFR"},
        ],
    )
    result = await fetch_pathways("P00533")
    assert len(result) == 2
    assert any("EGFR" in p["displayName"] for p in result)


@pytest.mark.asyncio
async def test_fetch_pathways_disease_only(httpx_mock: "HTTPXMock"):
    """disease_only=True should filter to disease pathways."""
    httpx_mock.add_response(
        url=f"{REACTOME_BASE}/data/pathways/low/entity/P13569?species=9606&disease=true",
        json=[{"stId": "R-HSA-5602410", "displayName": "CFTR in disease"}],
    )
    result = await fetch_pathways("P13569", disease_only=True)
    assert len(result) == 1
    assert "CFTR" in result[0]["displayName"]


@pytest.mark.asyncio
async def test_fetch_pathway_hierarchy(httpx_mock: "HTTPXMock"):
    """Should return ancestor pathway hierarchy as a list."""
    httpx_mock.add_response(
        url=f"{REACTOME_BASE}/data/pathway/R-HSA-177929/ancestors",
        json=[
            [{"stId": "R-HSA-162582", "displayName": "Signal Transduction"}],
            [{"stId": "R-HSA-1", "displayName": "Root"}],
        ],
    )
    result = await fetch_pathway_hierarchy("R-HSA-177929")
    assert isinstance(result, list)
    assert len(result) == 2


@pytest.mark.asyncio
async def test_fetch_contained_events(httpx_mock: "HTTPXMock"):
    """Should return reactions within a pathway."""
    httpx_mock.add_response(
        url=f"{REACTOME_BASE}/data/pathway/R-HSA-177929/containedEvents",
        json=[
            {"dbId": 12345, "displayName": "EGF binds EGFR"},
            {"dbId": 12346, "displayName": "EGFR phosphorylation"},
        ],
    )
    result = await fetch_contained_events("R-HSA-177929")
    assert len(result) == 2


@pytest.mark.asyncio
async def test_fetch_all_reactome_cftr(httpx_mock: "HTTPXMock"):
    """Integration: fetch_all_reactome should return both signaling and disease."""
    httpx_mock.add_response(
        url=f"{REACTOME_BASE}/data/pathways/low/entity/P13569?species=9606",
        json=[{"stId": "R-HSA-1234", "displayName": "CFTR transport"}],
    )
    httpx_mock.add_response(
        url=f"{REACTOME_BASE}/data/pathways/low/entity/P13569?species=9606&disease=true",
        json=[{"stId": "R-HSA-5602410", "displayName": "CFTR in disease"}],
    )
    result = await fetch_all_reactome("P13569")
    assert len(result["signaling_pathways"]) == 1
    assert len(result["disease_pathways"]) == 1


@pytest.mark.asyncio
async def test_fetch_all_reactome_bad_response(httpx_mock: "HTTPXMock"):
    """Non-200 responses should return empty lists, not raise exceptions."""
    httpx_mock.add_response(
        url=f"{REACTOME_BASE}/data/pathways/low/entity/INVALID?species=9606",
        status_code=404,
        text="Not found",
    )
    httpx_mock.add_response(
        url=f"{REACTOME_BASE}/data/pathways/low/entity/INVALID?species=9606&disease=true",
        status_code=404,
        text="Not found",
    )
    result = await fetch_all_reactome("INVALID")
    assert result["signaling_pathways"] == []
    assert result["disease_pathways"] == []
