"""Tests for analyzers/domain_extractor.py."""

from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from analyzers.domain_extractor import extract_features_from_uniprot, format_feature_summary
from core.models import SequenceFeature

CANONICAL = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"

FEATURES = [
    {
        "type": "Binding site",
        "location": {"start": {"value": 63}, "end": {"value": 63}},
        "description": "Heme iron coordination",
        "ligand": {"name": "Fe cation", "id": "ChEBI:CHEBI:18248"},
    },
    {
        "type": "Active site",
        "location": {"start": {"value": 92}, "end": {"value": 92}},
        "description": "Proton donor",
        "ligands": [],
    },
    {
        "type": "Domain",
        "location": {"start": {"value": 1}, "end": {"value": 147}},
        "description": "Globin",
        "ligands": [],
    },
    {
        "type": "Region",
        "location": {"start": {"value": 40}, "end": {"value": 50}},
        "description": "Interaction with alpha chain",
        "ligands": [],
    },
]


def test_extract_binding_sites():
    binding, active, domains = extract_features_from_uniprot(FEATURES, CANONICAL)
    assert len(binding) >= 1
    # The "Region" with "interaction" keyword should also be in binding
    types = [sf.feature_type for sf in binding]
    assert "Binding site" in types


def test_extract_active_sites():
    _, active, _ = extract_features_from_uniprot(FEATURES, CANONICAL)
    assert len(active) == 1
    assert active[0].start == 92


def test_extract_domains():
    _, _, domains = extract_features_from_uniprot(FEATURES, CANONICAL)
    domain_types = [sf.feature_type for sf in domains]
    assert "Domain" in domain_types


def test_sequence_fragment_slicing():
    binding, _, _ = extract_features_from_uniprot(FEATURES, CANONICAL)
    bs = next(sf for sf in binding if sf.feature_type == "Binding site")
    # Position 63 (1-based) => index 62
    assert bs.sequence_fragment == CANONICAL[62:63]


def test_ligand_extracted():
    binding, _, _ = extract_features_from_uniprot(FEATURES, CANONICAL)
    bs = next(sf for sf in binding if sf.feature_type == "Binding site")
    assert bs.ligand == "Fe cation"


def test_format_feature_summary():
    sf = SequenceFeature(
        feature_type="Binding site",
        start=63,
        end=63,
        sequence_fragment="H",
        description="Heme iron coordination",
        ligand="Fe cation",
    )
    summary = format_feature_summary(sf)
    assert "63" in summary
    assert "Fe cation" in summary
    assert "H" in summary
