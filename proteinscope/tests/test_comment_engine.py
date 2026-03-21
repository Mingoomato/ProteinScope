"""Tests for annotators/comment_engine.py"""

from __future__ import annotations

import pytest
from core.models import (
    DrugInteraction,
    DiseasePathway,
    MetabolicPathway,
    MetabolicReaction,
    PathologyStep,
    ProteinRecord,
    SignalingPathway,
)


def _minimal_record(**overrides) -> ProteinRecord:
    """Build a minimal ProteinRecord suitable for comment engine testing."""
    defaults = dict(
        query="EGFR",
        uniprot_id="P00533",
        protein_name="Epidermal growth factor receptor",
        gene_name="EGFR",
        organism="Homo sapiens",
        taxonomy_id="9606",
    )
    defaults.update(overrides)
    return ProteinRecord(**defaults)


# ---------------------------------------------------------------------------
# build_interaction_comment
# ---------------------------------------------------------------------------

def test_build_interaction_comment_metabolic():
    """Should use metabolic reaction data when available."""
    from annotators.comment_engine import build_interaction_comment

    record = _minimal_record(
        metabolic_pathways=[
            MetabolicPathway(
                pathway_id="hsa00010",
                pathway_name="Glycolysis",
                diagram_url="https://kegg.jp/pathway/hsa00010",
                reactions=[
                    MetabolicReaction(
                        reaction_id="R00001",
                        equation="ATP + Glucose -> ADP + Glucose-6P",
                        substrates=["ATP", "Glucose"],
                        products=["ADP", "Glucose-6P"],
                        cofactors_required=["Mg2+"],
                    )
                ],
            )
        ]
    )

    node = {"gene": "EGFR"}
    comment = build_interaction_comment(node, record)
    assert "EGFR" in comment
    assert "catalyzes" in comment
    assert "ATP" in comment or "Glucose" in comment


def test_build_interaction_comment_signaling():
    """Should fall back to signaling data when no metabolic reaction."""
    from annotators.comment_engine import build_interaction_comment

    record = _minimal_record(
        signaling_pathways=[
            SignalingPathway(
                pathway_id="R-HSA-177929",
                pathway_name="EGFR signaling",
                protein_role="receptor tyrosine kinase",
                activates=["RAS", "PI3K"],
                inhibits=["PTEN"],
                diagram_url="https://reactome.org/pathway/R-HSA-177929",
            )
        ]
    )

    node = {"gene": "EGFR"}
    comment = build_interaction_comment(node, record)
    assert "EGFR" in comment
    assert "receptor tyrosine kinase" in comment


def test_build_interaction_comment_fallback():
    """Should return generic fallback when no pathway data available."""
    from annotators.comment_engine import build_interaction_comment

    record = _minimal_record()
    node = {"gene": "EGFR"}
    comment = build_interaction_comment(node, record)
    assert "EGFR" in comment
    assert "participates" in comment


def test_build_interaction_comment_different_gene():
    """For a gene not matching record's gene_name, should still return something."""
    from annotators.comment_engine import build_interaction_comment

    record = _minimal_record()
    node = {"gene": "TP53"}
    comment = build_interaction_comment(node, record)
    assert "TP53" in comment


# ---------------------------------------------------------------------------
# build_drug_comment
# ---------------------------------------------------------------------------

def test_build_drug_comment_with_data():
    """Should list drug names when drug interactions exist for the gene."""
    from annotators.comment_engine import build_drug_comment

    record = _minimal_record(
        drug_interactions=[
            DrugInteraction(
                drug_name="Gefitinib",
                drug_id="CHEMBL939",
                interaction_type="inhibitor",
                clinical_significance="1A",
                sources=["ChEMBL"],
            ),
            DrugInteraction(
                drug_name="Erlotinib",
                drug_id="CHEMBL553",
                interaction_type="inhibitor",
                clinical_significance="1B",
                sources=["PharmGKB"],
            ),
        ]
    )

    node = {"gene": "EGFR"}
    comment = build_drug_comment(node, record)
    assert comment is not None
    assert "Gefitinib" in comment
    assert "Drug targets" in comment


def test_build_drug_comment_no_data():
    """Should return None when no drug interactions exist."""
    from annotators.comment_engine import build_drug_comment

    record = _minimal_record()
    node = {"gene": "EGFR"}
    comment = build_drug_comment(node, record)
    assert comment is None


# ---------------------------------------------------------------------------
# build_disease_comment
# ---------------------------------------------------------------------------

def test_build_disease_comment_with_data():
    """Should list disease names when disease pathway data exists."""
    from annotators.comment_engine import build_disease_comment

    record = _minimal_record(
        disease_pathways=[
            DiseasePathway(
                disease_name="Non-small cell lung cancer",
                pathway_id="R-HSA-5637815",
                pathway_source="Reactome",
                diagram_url="https://reactome.org",
                cascade=[
                    PathologyStep(
                        step_number=1,
                        event="EGFR mutation activates downstream signaling",
                        consequence="Uncontrolled proliferation",
                        therapeutic_target=True,
                    )
                ],
            )
        ]
    )

    node = {"gene": "EGFR"}
    comment = build_disease_comment(node, record)
    assert comment is not None
    assert "Non-small cell lung cancer" in comment
    assert "Disease association" in comment


def test_build_disease_comment_no_data():
    """Should return None when no disease pathway data."""
    from annotators.comment_engine import build_disease_comment

    record = _minimal_record()
    node = {"gene": "EGFR"}
    comment = build_disease_comment(node, record)
    assert comment is None


# ---------------------------------------------------------------------------
# build_stream_comment
# ---------------------------------------------------------------------------

def test_build_stream_comment_with_data():
    """Should return upstream/downstream strings when signaling data exists."""
    from annotators.comment_engine import build_stream_comment

    record = _minimal_record(
        signaling_pathways=[
            SignalingPathway(
                pathway_id="R-HSA-177929",
                pathway_name="EGFR signaling",
                activated_by=["EGF", "TGF-alpha"],
                activates=["RAS", "PI3K", "STAT3"],
                diagram_url="https://reactome.org",
            )
        ]
    )

    node = {"gene": "EGFR"}
    upstream, downstream = build_stream_comment(node, record)
    assert upstream is not None
    assert "EGF" in upstream
    assert downstream is not None
    assert "RAS" in downstream


def test_build_stream_comment_no_data():
    """Should return (None, None) when no signaling data."""
    from annotators.comment_engine import build_stream_comment

    record = _minimal_record()
    node = {"gene": "EGFR"}
    upstream, downstream = build_stream_comment(node, record)
    assert upstream is None
    assert downstream is None


# ---------------------------------------------------------------------------
# build_arrow_label
# ---------------------------------------------------------------------------

def test_build_arrow_label_with_drug():
    """Should prioritize drug info in arrow label."""
    from annotators.comment_engine import build_arrow_label, MAX_ARROW_CHARS

    interaction = "EGFR catalyzes phosphorylation"
    drug = "Drug targets: Gefitinib (inhibitor, 1A)"
    label = build_arrow_label(interaction, drug)
    assert "Drug target" in label
    assert "Gefitinib" in label
    assert len(label) <= MAX_ARROW_CHARS


def test_build_arrow_label_without_drug():
    """Should use first words of interaction comment when no drug data."""
    from annotators.comment_engine import build_arrow_label, MAX_ARROW_CHARS

    interaction = "EGFR acts as receptor tyrosine kinase and phosphorylates downstream targets"
    label = build_arrow_label(interaction, None)
    assert len(label) <= MAX_ARROW_CHARS
    assert "EGFR" in label


def test_build_arrow_label_truncation():
    """Long labels should be truncated to MAX_ARROW_CHARS."""
    from annotators.comment_engine import build_arrow_label, MAX_ARROW_CHARS

    long_interaction = " ".join(["word"] * 50)
    label = build_arrow_label(long_interaction, None)
    assert len(label) <= MAX_ARROW_CHARS


# ---------------------------------------------------------------------------
# get_protein_name
# ---------------------------------------------------------------------------

def test_get_protein_name_same_gene():
    """Should return record's protein_name when gene matches."""
    from annotators.comment_engine import get_protein_name

    record = _minimal_record()
    name = get_protein_name("EGFR", record)
    assert name == "Epidermal growth factor receptor"


def test_get_protein_name_case_insensitive():
    """Gene matching should be case-insensitive."""
    from annotators.comment_engine import get_protein_name

    record = _minimal_record()
    name = get_protein_name("egfr", record)
    assert name == "Epidermal growth factor receptor"


def test_get_protein_name_different_gene():
    """Should return the gene symbol itself when gene doesn't match record."""
    from annotators.comment_engine import get_protein_name

    record = _minimal_record()
    name = get_protein_name("TP53", record)
    assert name == "TP53"
