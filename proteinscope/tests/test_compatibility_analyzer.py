"""Tests for compatibility_analyzer.py"""
import pytest
from analyzers.compatibility_analyzer import CompatibilityReport, run_compatibility_analysis


def test_compatibility_report_model():
    report = CompatibilityReport(
        source_gene="POLD1",
        source_organism="penguin",
        target_gene="POLD1",
        target_organism="Homo sapiens",
        target_cell="HeLa",
    )
    assert report.source_gene == "POLD1"
    assert report.verdict == ""
    assert report.sequence_identity_pct is None


def test_compatibility_report_with_data():
    report = CompatibilityReport(
        source_gene="POLD1",
        source_organism="penguin",
        target_gene="POLD1",
        target_organism="Homo sapiens",
        sequence_identity_pct=91.2,
        tm_score=0.94,
        verdict="Likely compatible",
        confidence="Moderate",
        reasoning="High sequence identity suggests functional conservation.",
    )
    assert report.sequence_identity_pct == 91.2
    assert report.verdict == "Likely compatible"


def test_compatibility_report_defaults():
    report = CompatibilityReport(
        source_gene="EGFR",
        source_organism="mouse",
        target_gene="EGFR",
        target_organism="Homo sapiens",
    )
    assert report.target_cell is None
    assert report.caveats == []
    assert report.recommendations == []
    assert report.source_uniprot_id is None
    assert report.target_uniprot_id is None
    assert report.confidence == ""


def test_compatibility_report_caveats_and_recommendations():
    report = CompatibilityReport(
        source_gene="TP53",
        source_organism="mouse",
        target_gene="TP53",
        target_organism="Homo sapiens",
        verdict="Conditional",
        confidence="Moderate",
        caveats=["Active site residue Q144 differs", "Temperature optimum 2°C lower"],
        recommendations=["Validate in cell-free transcription assay", "Check p53 response element binding"],
    )
    assert len(report.caveats) == 2
    assert len(report.recommendations) == 2
    assert "Q144" in report.caveats[0]


def test_compatibility_report_active_site_conservation():
    report = CompatibilityReport(
        source_gene="PCNA",
        source_organism="zebrafish",
        target_gene="PCNA",
        target_organism="Homo sapiens",
        active_site_conservation_pct=100.0,
        sequence_identity_pct=96.5,
        verdict="Compatible",
        confidence="High",
    )
    assert report.active_site_conservation_pct == 100.0
    assert report.verdict == "Compatible"
