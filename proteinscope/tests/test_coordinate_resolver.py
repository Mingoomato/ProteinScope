"""Tests for annotators/coordinate_resolver.py"""

from __future__ import annotations

import asyncio
from unittest.mock import AsyncMock, MagicMock, patch

import pytest


# ---------------------------------------------------------------------------
# KEGG coordinate tests
# ---------------------------------------------------------------------------

KGML_SAMPLE = b"""<?xml version="1.0"?>
<pathway name="path:hsa00010" number="00010" title="Glycolysis">
  <image name="hsa00010" width="1024" height="768"/>
  <entry id="1" type="gene">
    <graphics name="EGFR, 1956" x="200" y="150" width="46" height="17"/>
  </entry>
  <entry id="2" type="gene">
    <graphics name="TP53, 7157" x="400" y="300" width="46" height="17"/>
  </entry>
</pathway>
"""


@pytest.mark.asyncio
async def test_get_kegg_node_coordinates_found():
    """Should return coordinate dict for a gene found in KGML."""
    from annotators.coordinate_resolver import get_kegg_node_coordinates

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.content = KGML_SAMPLE

    with patch("httpx.AsyncClient") as mock_client_cls:
        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.get = AsyncMock(return_value=mock_response)
        mock_client_cls.return_value = mock_client

        results = await get_kegg_node_coordinates(
            "hsa00010", ["EGFR"], image_width=1024, image_height=768
        )

    assert len(results) == 1
    r = results[0]
    assert r["gene"] == "EGFR"
    assert r["source"] == "KEGG"
    assert r["pathway_id"] == "hsa00010"
    assert r["x"] == 200
    assert r["y"] == 150


@pytest.mark.asyncio
async def test_get_kegg_node_coordinates_not_found():
    """Should return empty list when gene not in KGML."""
    from annotators.coordinate_resolver import get_kegg_node_coordinates

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.content = KGML_SAMPLE

    with patch("httpx.AsyncClient") as mock_client_cls:
        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.get = AsyncMock(return_value=mock_response)
        mock_client_cls.return_value = mock_client

        results = await get_kegg_node_coordinates(
            "hsa00010", ["BRCA1"], image_width=1024, image_height=768
        )

    assert results == []


@pytest.mark.asyncio
async def test_get_kegg_node_coordinates_http_error():
    """Should return empty list on non-200 response."""
    from annotators.coordinate_resolver import get_kegg_node_coordinates

    mock_response = MagicMock()
    mock_response.status_code = 404

    with patch("httpx.AsyncClient") as mock_client_cls:
        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.get = AsyncMock(return_value=mock_response)
        mock_client_cls.return_value = mock_client

        results = await get_kegg_node_coordinates("hsa99999", ["EGFR"])

    assert results == []


@pytest.mark.asyncio
async def test_get_kegg_node_coordinates_scale():
    """Scale should be applied when image dimensions differ from KGML canvas."""
    from annotators.coordinate_resolver import get_kegg_node_coordinates

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.content = KGML_SAMPLE  # KGML canvas: 1024x768

    with patch("httpx.AsyncClient") as mock_client_cls:
        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.get = AsyncMock(return_value=mock_response)
        mock_client_cls.return_value = mock_client

        # Request with double the canvas size → coordinates should be doubled
        results = await get_kegg_node_coordinates(
            "hsa00010", ["EGFR"], image_width=2048, image_height=1536
        )

    assert len(results) == 1
    r = results[0]
    assert r["x"] == 400   # 200 * (2048/1024)
    assert r["y"] == 300   # 150 * (1536/768)


@pytest.mark.asyncio
async def test_get_kegg_node_coordinates_network_exception():
    """Should return empty list on network exception."""
    from annotators.coordinate_resolver import get_kegg_node_coordinates

    with patch("httpx.AsyncClient") as mock_client_cls:
        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.get = AsyncMock(side_effect=Exception("timeout"))
        mock_client_cls.return_value = mock_client

        results = await get_kegg_node_coordinates("hsa00010", ["EGFR"])

    assert results == []


# ---------------------------------------------------------------------------
# Reactome coordinate tests
# ---------------------------------------------------------------------------

REACTOME_LAYOUT = {
    "minX": 0, "minY": 0, "maxX": 800, "maxY": 600,
    "nodes": [
        {
            "displayName": "EGFR protein",
            "schemaClass": "Protein",
            "stId": "R-HSA-1234",
            "x": 200.0,
            "y": 150.0,
            "width": 60,
            "height": 25,
        },
        {
            "displayName": "Some complex",
            "schemaClass": "Complex",
            "stId": "R-HSA-5678",
            "x": 400.0,
            "y": 300.0,
            "width": 80,
            "height": 30,
        },
    ]
}


@pytest.mark.asyncio
async def test_get_reactome_node_coordinates_found():
    """Should return coordinate dict for a gene found in Reactome layout."""
    from annotators.coordinate_resolver import get_reactome_node_coordinates

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json = MagicMock(return_value=REACTOME_LAYOUT)

    with patch("httpx.AsyncClient") as mock_client_cls:
        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.get = AsyncMock(return_value=mock_response)
        mock_client_cls.return_value = mock_client

        results = await get_reactome_node_coordinates(
            "R-HSA-177929", ["EGFR"], image_width=800, image_height=600
        )

    assert len(results) == 1
    r = results[0]
    assert r["gene"] == "EGFR"
    assert r["source"] == "Reactome"
    assert r["x"] == 200
    assert r["y"] == 150


@pytest.mark.asyncio
async def test_get_reactome_node_coordinates_schema_filter():
    """Should only match nodes with valid schemaClass values."""
    from annotators.coordinate_resolver import get_reactome_node_coordinates

    layout = {
        "minX": 0, "minY": 0, "maxX": 800, "maxY": 600,
        "nodes": [
            {
                "displayName": "EGFR mRNA",
                "schemaClass": "RNA",  # should be filtered out
                "stId": "R-HSA-999",
                "x": 100.0, "y": 100.0, "width": 50, "height": 20,
            }
        ]
    }

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json = MagicMock(return_value=layout)

    with patch("httpx.AsyncClient") as mock_client_cls:
        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.get = AsyncMock(return_value=mock_response)
        mock_client_cls.return_value = mock_client

        results = await get_reactome_node_coordinates(
            "R-HSA-177929", ["EGFR"], image_width=800, image_height=600
        )

    assert results == []


# ---------------------------------------------------------------------------
# PharmGKB fallback test
# ---------------------------------------------------------------------------

@pytest.mark.asyncio
async def test_pharmgkb_fallback_on_404():
    """Should return fallback grid positions when SVG download fails."""
    from annotators.coordinate_resolver import get_pharmgkb_node_coordinates

    mock_response = MagicMock()
    mock_response.status_code = 404

    with patch("httpx.AsyncClient") as mock_client_cls:
        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.get = AsyncMock(return_value=mock_response)
        mock_client_cls.return_value = mock_client

        results = await get_pharmgkb_node_coordinates(
            "PA150654557", ["CYP2C9", "VKORC1"]
        )

    assert len(results) == 2
    for r in results:
        assert r["source"] == "PharmGKB"
        assert r["x"] >= 10
        assert r["y"] >= 10


# ---------------------------------------------------------------------------
# STRING coordinate tests
# ---------------------------------------------------------------------------

@pytest.mark.asyncio
async def test_get_string_node_coordinates_found():
    """Should return position for the queried gene from STRING API."""
    from annotators.coordinate_resolver import get_string_node_coordinates

    string_data = [
        {
            "preferredName": "EGFR",
            "stringId": "9606.ENSP00000275493",
            "x": 512.0,
            "y": 384.0,
        }
    ]

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json = MagicMock(return_value=string_data)

    with patch("httpx.AsyncClient") as mock_client_cls:
        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.get = AsyncMock(return_value=mock_response)
        mock_client_cls.return_value = mock_client

        results = await get_string_node_coordinates("EGFR")

    assert len(results) == 1
    r = results[0]
    assert r["gene"] == "EGFR"
    assert r["source"] == "STRING"
    assert r["x"] == 512
    assert r["y"] == 384


@pytest.mark.asyncio
async def test_get_string_node_coordinates_empty():
    """Should return empty list when API returns non-list."""
    from annotators.coordinate_resolver import get_string_node_coordinates

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json = MagicMock(return_value={"error": "not found"})

    with patch("httpx.AsyncClient") as mock_client_cls:
        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.get = AsyncMock(return_value=mock_response)
        mock_client_cls.return_value = mock_client

        results = await get_string_node_coordinates("UNKNOWNGENE")

    assert results == []
