"""Tests for fetchers/pharmgkb.py using pytest-httpx."""

from __future__ import annotations

import json
import pytest

try:
    from pytest_httpx import HTTPXMock
    HAS_HTTPX_MOCK = True
except ImportError:
    HAS_HTTPX_MOCK = False

from fetchers.pharmgkb import (
    get_pharmgkb_id,
    fetch_drug_interactions,
    fetch_pgx_variants,
    fetch_pharmgkb_pathways,
    fetch_fda_labels,
    fetch_all_pharmgkb,
)

pytestmark = pytest.mark.skipif(
    not HAS_HTTPX_MOCK,
    reason="pytest-httpx not installed",
)

BASE = "https://api.pharmgkb.org/v1/data"


@pytest.mark.asyncio
async def test_get_pharmgkb_id_cyp2c9(httpx_mock: "HTTPXMock"):
    """CYP2C9 should resolve to PharmGKB gene ID PA126."""
    httpx_mock.add_response(
        url=f"{BASE}/gene?symbol=CYP2C9&view=min",
        json={"data": [{"id": "PA126", "symbol": "CYP2C9"}]},
    )
    result = await get_pharmgkb_id("CYP2C9")
    assert result == "PA126"


@pytest.mark.asyncio
async def test_get_pharmgkb_id_not_found(httpx_mock: "HTTPXMock"):
    """Unknown gene returns None."""
    httpx_mock.add_response(
        url=f"{BASE}/gene?symbol=NOTREAL&view=min",
        json={"data": []},
    )
    result = await get_pharmgkb_id("NOTREAL")
    assert result is None


@pytest.mark.asyncio
async def test_fetch_drug_interactions_cyp2c9(httpx_mock: "HTTPXMock"):
    """Drug interactions for CYP2C9 (PA126) should include Warfarin."""
    httpx_mock.add_response(
        url=f"{BASE}/gene/PA126/relatedChemicals",
        json={"data": [{"id": "PA451906", "name": "warfarin", "type": "substrate"}]},
    )
    result = await fetch_drug_interactions("PA126")
    assert len(result) == 1
    assert result[0]["name"] == "warfarin"


@pytest.mark.asyncio
async def test_fetch_pgx_variants_cyp2c9(httpx_mock: "HTTPXMock"):
    """PGx variants for CYP2C9 should include rs1799853."""
    httpx_mock.add_response(
        url=f"{BASE}/clinicalAnnotation?gene=PA126&view=max",
        json={
            "data": [
                {
                    "evidenceLevel": "1A",
                    "location": {"variants": [{"name": "rs1799853"}]},
                    "relatedChemicals": [{"name": "warfarin"}],
                    "phenotype": "reduced metabolism",
                }
            ]
        },
    )
    result = await fetch_pgx_variants("PA126")
    assert len(result) == 1
    assert result[0]["evidenceLevel"] == "1A"


@pytest.mark.asyncio
async def test_fetch_pharmgkb_pathways(httpx_mock: "HTTPXMock"):
    """Should return pathway list for CYP2C9."""
    httpx_mock.add_response(
        url=f"{BASE}/pathway?relatedGene=PA126",
        json={"data": [{"id": "PA150654557", "name": "Warfarin Pathway, Pharmacodynamics"}]},
    )
    result = await fetch_pharmgkb_pathways("PA126")
    assert len(result) == 1
    assert "Warfarin" in result[0]["name"]


@pytest.mark.asyncio
async def test_fetch_fda_labels_cyp2c9(httpx_mock: "HTTPXMock"):
    """Should return FDA labels mentioning CYP2C9."""
    httpx_mock.add_response(
        url=f"{BASE}/label?relatedGene=PA126",
        json={"data": [{"id": "L01", "name": "Warfarin FDA Label"}]},
    )
    result = await fetch_fda_labels("PA126")
    assert len(result) == 1


@pytest.mark.asyncio
async def test_fetch_all_pharmgkb_cyp2c9(httpx_mock: "HTTPXMock"):
    """Integration: fetch_all_pharmgkb should return structured dict."""
    httpx_mock.add_response(
        url=f"{BASE}/gene?symbol=CYP2C9&view=min",
        json={"data": [{"id": "PA126", "symbol": "CYP2C9"}]},
    )
    httpx_mock.add_response(
        url=f"{BASE}/gene/PA126/relatedChemicals",
        json={"data": [{"id": "PA451906", "name": "warfarin"}]},
    )
    httpx_mock.add_response(
        url=f"{BASE}/clinicalAnnotation?gene=PA126&view=max",
        json={"data": []},
    )
    httpx_mock.add_response(
        url=f"{BASE}/pathway?relatedGene=PA126",
        json={"data": []},
    )
    httpx_mock.add_response(
        url=f"{BASE}/label?relatedGene=PA126",
        json={"data": []},
    )
    result = await fetch_all_pharmgkb("CYP2C9")
    assert result["pharmgkb_id"] == "PA126"
    assert len(result["drug_interactions"]) == 1
