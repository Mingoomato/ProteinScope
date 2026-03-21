"""Tests for fetchers/uniprot.py using pytest-httpx fixtures."""

from __future__ import annotations

import json
import pytest
from pytest_httpx import HTTPXMock

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from fetchers.uniprot import (
    fetch_by_accession,
    extract_protein_name,
    extract_gene_name,
    extract_organism,
    extract_function,
    extract_cofactors,
    extract_subcellular_location,
    extract_features,
    extract_references,
)

# ── Fixtures ──────────────────────────────────────────────────────────────────

HEMOGLOBIN_ENTRY = {
    "primaryAccession": "P68871",
    "proteinDescription": {
        "recommendedName": {
            "fullName": {"value": "Hemoglobin subunit beta"},
            "ecNumbers": [],
        }
    },
    "genes": [{"geneName": {"value": "HBB"}}],
    "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
    "sequence": {"value": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH", "length": 147},
    "comments": [
        {
            "commentType": "FUNCTION",
            "texts": [{"value": "Involved in oxygen transport from the lung to the various peripheral tissues."}],
        },
        {
            "commentType": "COFACTOR",
            "cofactors": [{"name": "Heme"}],
        },
        {
            "commentType": "SUBCELLULAR LOCATION",
            "subcellularLocations": [{"location": {"value": "Cytoplasm"}}],
        },
    ],
    "features": [
        {
            "type": "Binding site",
            "location": {"start": {"value": 63}, "end": {"value": 63}},
            "description": "Heme iron coordination",
            "ligands": [{"name": "Fe cation"}],
        }
    ],
    "references": [
        {
            "citation": {
                "title": "Sequence of amino acids in the alpha and beta chains of adult human hemoglobin.",
                "authors": [{"value": "Ingram VM"}],
                "journal": "Nature",
                "publicationDate": "1956",
                "citationCrossReferences": [
                    {"database": "PubMed", "id": "13369412"},
                ],
            }
        }
    ],
    "uniProtKBCrossReferences": [],
}


@pytest.mark.asyncio
async def test_fetch_by_accession(httpx_mock: HTTPXMock):
    httpx_mock.add_response(
        url="https://rest.uniprot.org/uniprotkb/P68871?format=json",
        json=HEMOGLOBIN_ENTRY,
    )
    result = await fetch_by_accession("P68871")
    assert result["primaryAccession"] == "P68871"


def test_extract_protein_name():
    assert extract_protein_name(HEMOGLOBIN_ENTRY) == "Hemoglobin subunit beta"


def test_extract_gene_name():
    assert extract_gene_name(HEMOGLOBIN_ENTRY) == "HBB"


def test_extract_organism():
    name, taxid = extract_organism(HEMOGLOBIN_ENTRY)
    assert name == "Homo sapiens"
    assert taxid == "9606"


def test_extract_function():
    fn = extract_function(HEMOGLOBIN_ENTRY)
    assert "oxygen transport" in fn.lower()


def test_extract_cofactors():
    cofs = extract_cofactors(HEMOGLOBIN_ENTRY)
    assert "Heme" in cofs


def test_extract_subcellular_location():
    locs = extract_subcellular_location(HEMOGLOBIN_ENTRY)
    assert "Cytoplasm" in locs


def test_extract_features():
    seq = HEMOGLOBIN_ENTRY["sequence"]["value"]
    binding, active, domains = extract_features(HEMOGLOBIN_ENTRY, seq)
    assert len(binding) == 1
    assert binding[0]["start"] == 63
    assert binding[0]["ligand"] == "Fe cation"


def test_extract_references():
    refs = extract_references(HEMOGLOBIN_ENTRY)
    assert len(refs) == 1
    assert refs[0]["pubmed_id"] == "13369412"
    assert "Ingram" in refs[0]["authors"][0]
