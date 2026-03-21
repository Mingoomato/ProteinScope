"""Tests for analyzers/pathway_integrator.py."""

from __future__ import annotations

import pytest

from analyzers.pathway_integrator import (
    classify_pathway,
    build_metabolic_pathways_from_kegg,
    build_signaling_pathways_from_reactome,
    build_disease_pathways_from_reactome,
    build_pharmacogenomic_pathways_from_pharmgkb,
    build_protein_interactions_from_string,
)
from analyzers.pgx_reporter import build_drug_interactions, build_pgx_variants


# ---------------------------------------------------------------------------
# classify_pathway
# ---------------------------------------------------------------------------

class TestClassifyPathway:
    def test_glycolysis(self):
        assert classify_pathway("Glycolysis / Gluconeogenesis") == "metabolism"

    def test_tca(self):
        assert classify_pathway("Citric acid cycle (TCA cycle)") == "metabolism"

    def test_mapk(self):
        assert classify_pathway("MAPK signaling pathway") == "signaling"

    def test_pi3k(self):
        assert classify_pathway("PI3K-Akt signaling pathway") == "signaling"

    def test_cancer(self):
        assert classify_pathway("Pathways in cancer") == "disease"

    def test_cyp450(self):
        assert classify_pathway("Drug metabolism - cytochrome P450") == "pharmacology"

    def test_other(self):
        assert classify_pathway("Completely unknown process") == "other"


# ---------------------------------------------------------------------------
# build_metabolic_pathways_from_kegg
# ---------------------------------------------------------------------------

class TestBuildMetabolicPathways:
    def test_glycolysis_classified(self):
        kegg_data = {
            "pathway_ids": ["path:hsa00010"],
            "pathway_details": {
                "path:hsa00010": {"NAME": "Glycolysis / Gluconeogenesis", "ENTRY": "hsa00010"}
            },
        }
        pathways = build_metabolic_pathways_from_kegg(kegg_data, "LDHA")
        assert len(pathways) == 1
        assert pathways[0].pathway_id == "hsa00010"
        assert "Glycolysis" in pathways[0].pathway_name

    def test_signaling_pathway_excluded(self):
        kegg_data = {
            "pathway_ids": ["path:hsa04010"],
            "pathway_details": {
                "path:hsa04010": {"NAME": "MAPK signaling pathway", "ENTRY": "hsa04010"}
            },
        }
        # MAPK is signaling, not metabolism — should be excluded from metabolic results
        pathways = build_metabolic_pathways_from_kegg(kegg_data, "EGFR")
        assert len(pathways) == 0

    def test_diagram_url_format(self):
        kegg_data = {
            "pathway_ids": ["path:hsa00020"],
            "pathway_details": {
                "path:hsa00020": {"NAME": "Citric acid cycle (TCA cycle)", "ENTRY": "hsa00020"}
            },
        }
        pathways = build_metabolic_pathways_from_kegg(kegg_data, "CS")
        assert "genome.jp" in pathways[0].diagram_url


# ---------------------------------------------------------------------------
# build_signaling_pathways_from_reactome
# ---------------------------------------------------------------------------

class TestBuildSignalingPathways:
    def test_egfr_signaling(self):
        reactome_data = {
            "signaling_pathways": [
                {"stId": "R-HSA-177929", "displayName": "Signaling by EGFR"},
            ],
            "disease_pathways": [],
        }
        pathways = build_signaling_pathways_from_reactome(reactome_data, "P00533", "EGFR")
        assert len(pathways) == 1
        assert pathways[0].pathway_id == "R-HSA-177929"

    def test_deduplication(self):
        reactome_data = {
            "signaling_pathways": [
                {"stId": "R-HSA-177929", "displayName": "Signaling by EGFR"},
                {"stId": "R-HSA-177929", "displayName": "Signaling by EGFR"},  # duplicate
            ],
            "disease_pathways": [],
        }
        pathways = build_signaling_pathways_from_reactome(reactome_data, "P00533", "EGFR")
        assert len(pathways) == 1

    def test_disease_pathways_excluded(self):
        reactome_data = {
            "signaling_pathways": [
                {"stId": "R-HSA-5602410", "displayName": "CFTR in cancer"}
            ],
            "disease_pathways": [],
        }
        # "cancer" in name → classified as disease → excluded from signaling
        pathways = build_signaling_pathways_from_reactome(reactome_data, "P13569", "CFTR")
        assert len(pathways) == 0


# ---------------------------------------------------------------------------
# build_disease_pathways_from_reactome
# ---------------------------------------------------------------------------

class TestBuildDiseasePathways:
    def test_cftr_disease(self):
        reactome_data = {
            "signaling_pathways": [],
            "disease_pathways": [
                {"stId": "R-HSA-5602410", "displayName": "CFTR in disease"}
            ],
        }
        pathways = build_disease_pathways_from_reactome(reactome_data, {}, "CFTR")
        assert len(pathways) == 1
        assert pathways[0].pathway_source == "Reactome"
        assert "reactome.org" in pathways[0].diagram_url

    def test_kegg_disease_ids_appended(self):
        reactome_data = {"signaling_pathways": [], "disease_pathways": []}
        kegg_data = {"disease_pathway_ids": ["ds:H00218"]}
        pathways = build_disease_pathways_from_reactome(reactome_data, kegg_data, "CFTR")
        assert len(pathways) == 1
        assert pathways[0].pathway_source == "KEGG"


# ---------------------------------------------------------------------------
# build_pharmacogenomic_pathways_from_pharmgkb
# ---------------------------------------------------------------------------

class TestBuildPharmaGenomicPathways:
    def test_warfarin_pathway(self):
        pharmgkb_data = {
            "pathways": [
                {
                    "id": "PA150654557",
                    "name": "Warfarin Pathway, Pharmacodynamics",
                    "relatedChemicals": [{"name": "warfarin"}],
                    "relatedGenes": [{"symbol": "CYP2C9"}],
                }
            ]
        }
        pathways = build_pharmacogenomic_pathways_from_pharmgkb(pharmgkb_data)
        assert len(pathways) == 1
        assert pathways[0].pathway_id == "PA150654557"
        assert "warfarin" in pathways[0].drugs_involved
        assert "CYP2C9" in pathways[0].genes_involved


# ---------------------------------------------------------------------------
# build_protein_interactions_from_string
# ---------------------------------------------------------------------------

class TestBuildProteinInteractions:
    def test_score_conversion(self):
        string_data = {
            "interactions": [
                {"preferredName_B": "TP53", "score": 900, "escore": 100},
                {"preferredName_B": "MDM2", "score": 300, "escore": 0},  # below default threshold
            ]
        }
        interactions = build_protein_interactions_from_string(string_data, min_score=0.4)
        assert len(interactions) == 1
        assert interactions[0].partner_protein == "TP53"
        assert interactions[0].confidence_score == pytest.approx(0.9)
        assert interactions[0].experimentally_confirmed is True

    def test_experimental_flag(self):
        string_data = {
            "interactions": [
                {"preferredName_B": "KRAS", "score": 600, "escore": 0},
            ]
        }
        interactions = build_protein_interactions_from_string(string_data)
        assert interactions[0].experimentally_confirmed is False


# ---------------------------------------------------------------------------
# build_drug_interactions (pgx_reporter)
# ---------------------------------------------------------------------------

class TestBuildDrugInteractions:
    def test_merges_sources(self):
        pharmgkb_data = {
            "drug_interactions": [{"id": "PA1", "name": "warfarin", "type": "substrate"}],
            "pgx_variants": [],
            "pathways": [],
            "fda_labels": [],
        }
        dgidb = [
            {"drug": {"name": "warfarin", "conceptId": "DRUGBANK:DB00682"},
             "interactionTypes": [{"type": "substrate", "directionality": None}],
             "interactionScore": 5.0, "publications": [], "sources": []}
        ]
        chembl_data = {
            "mechanisms": [
                {"molecule_name": "warfarin", "molecule_chembl_id": "CHEMBL1599",
                 "action_type": "SUBSTRATE", "mechanism_of_action": "CYP2C9 hydroxylation"}
            ],
            "activities": [],
        }
        result = build_drug_interactions(pharmgkb_data, dgidb, chembl_data)
        assert len(result) == 1
        di = result[0]
        assert "PharmGKB" in di.sources
        assert "DGIdb" in di.sources
        assert "ChEMBL" in di.sources
        assert di.mechanism == "CYP2C9 hydroxylation"

    def test_fda_label_flag(self):
        pharmgkb_data = {
            "drug_interactions": [{"id": "PA1", "name": "warfarin", "type": "substrate"}],
            "pgx_variants": [],
            "pathways": [],
            "fda_labels": [{"relatedChemicals": [{"name": "warfarin"}]}],
        }
        result = build_drug_interactions(pharmgkb_data, [], {})
        assert result[0].fda_label_mentioned is True


# ---------------------------------------------------------------------------
# build_pgx_variants (pgx_reporter)
# ---------------------------------------------------------------------------

class TestBuildPGxVariants:
    def test_evidence_level_filter(self):
        pharmgkb_data = {
            "pgx_variants": [
                {
                    "evidenceLevel": "1A",
                    "location": {"variants": [{"name": "rs1799853"}]},
                    "relatedChemicals": [{"name": "warfarin"}],
                    "phenotype": "reduced metabolism",
                    "evidence": [],
                },
                {
                    "evidenceLevel": "4",
                    "location": {"variants": [{"name": "rs9999"}]},
                    "relatedChemicals": [{"name": "aspirin"}],
                    "phenotype": "mild effect",
                    "evidence": [],
                },
            ]
        }
        # min_level 1B → only 1A passes
        result = build_pgx_variants(pharmgkb_data, "CYP2C9", min_pgx_evidence="1B")
        assert len(result) == 1
        assert result[0].variant_id == "rs1799853"

    def test_all_evidence_levels_included_by_default(self):
        pharmgkb_data = {
            "pgx_variants": [
                {
                    "evidenceLevel": "4",
                    "location": {"variants": [{"name": "rs123"}]},
                    "relatedChemicals": [{"name": "aspirin"}],
                    "phenotype": "effect",
                    "evidence": [],
                }
            ]
        }
        result = build_pgx_variants(pharmgkb_data, "GENE1", min_pgx_evidence="4")
        assert len(result) == 1
