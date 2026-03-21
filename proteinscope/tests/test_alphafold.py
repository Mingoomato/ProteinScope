"""Tests for fetchers/alphafold.py."""

from __future__ import annotations

import pytest
from pytest_httpx import HTTPXMock

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from fetchers.alphafold import fetch_alphafold_metadata

MOCK_METADATA = [
    {
        "entryId": "AF-P68871-F1",
        "gene": "HBB",
        "uniprotAccession": "P68871",
        "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P68871-F1-model_v4.pdb",
        "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P68871-F1-model_v4.cif",
        "confidenceScore": 92.3,
        "organismScientificName": "Homo sapiens",
    }
]


@pytest.mark.asyncio
async def test_fetch_alphafold_metadata(httpx_mock: HTTPXMock):
    httpx_mock.add_response(
        url="https://alphafold.ebi.ac.uk/api/prediction/P68871",
        json=MOCK_METADATA,
    )
    meta = await fetch_alphafold_metadata("P68871")
    assert meta is not None
    assert meta["uniprotAccession"] == "P68871"
    assert meta["confidenceScore"] == 92.3
    assert "pdbUrl" in meta


@pytest.mark.asyncio
async def test_fetch_alphafold_metadata_not_found(httpx_mock: HTTPXMock):
    httpx_mock.add_response(
        url="https://alphafold.ebi.ac.uk/api/prediction/INVALID",
        status_code=404,
    )
    meta = await fetch_alphafold_metadata("INVALID")
    assert meta is None
