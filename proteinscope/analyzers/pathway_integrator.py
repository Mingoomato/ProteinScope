"""Pathway integrator — merge KEGG + Reactome + MetaCyc results.

Deduplicates overlapping pathway entries and classifies each into:
  metabolism | signaling | disease | pharmacology | other
"""

from __future__ import annotations

from core.models import (
    DiseasePathway,
    MetabolicPathway,
    MetabolicReaction,
    PathologyStep,
    PharmacogenomicPathway,
    ProteinInteraction,
    SignalingPathway,
)

# ---------------------------------------------------------------------------
# Pathway classification
# ---------------------------------------------------------------------------

_CATEGORY_MAP = {
    "metabolism": [
        "glycolysis", "tca", "citric acid", "krebs", "fatty acid",
        "amino acid", "nucleotide", "oxidative phosphorylation",
        "pentose phosphate", "gluconeogenesis", "lipid", "steroid",
        "pyruvate", "carbon", "nitrogen",
    ],
    "signaling": [
        "mapk", "pi3k", "akt", "wnt", "notch", "hedgehog",
        "jak-stat", "mtor", "camp", "calcium", "pkc",
        "nf-kb", "tgf", "egfr", "insulin signaling",
        "apoptosis", "cell cycle", "p53", "ras",
    ],
    "disease": [
        "cancer", "tumor", "carcinoma", "leukemia",
        "diabetes", "obesity", "neurodegeneration", "alzheimer",
        "parkinson", "huntington", "infection", "immune disorder",
        "cystic fibrosis", "asthma", "hypertension",
    ],
    "pharmacology": [
        "drug metabolism", "cyp450", "cytochrome p450", "transport",
        "pharmacokinetics", "abc transporter", "solute carrier",
        "xenobiotic", "bile acid",
    ],
}


def classify_pathway(pathway_name: str) -> str:
    """Return the primary category for a pathway based on its name."""
    name_lower = pathway_name.lower()
    for category, keywords in _CATEGORY_MAP.items():
        if any(kw in name_lower for kw in keywords):
            return category
    return "other"


# ---------------------------------------------------------------------------
# KEGG pathway assembly
# ---------------------------------------------------------------------------

def build_metabolic_pathways_from_kegg(kegg_data: dict, gene_name: str) -> list[MetabolicPathway]:
    """Convert KEGG fetch results into MetabolicPathway models."""
    pathways = []
    for pid, detail in kegg_data.get("pathway_details", {}).items():
        name = detail.get("NAME", "").strip() or pid
        category = classify_pathway(name)
        if category not in ("metabolism", "pharmacology", "other"):
            continue

        clean_pid = pid.replace("path:", "")
        pathways.append(MetabolicPathway(
            pathway_id=clean_pid,
            pathway_name=name,
            protein_role=gene_name,
            reactions=[],
            upstream_proteins=[],
            downstream_proteins=[],
            diagram_url=f"https://www.genome.jp/pathway/{clean_pid}",
            diagram_image_path=None,
        ))
    return pathways


# ---------------------------------------------------------------------------
# Reactome pathway assembly
# ---------------------------------------------------------------------------

def build_signaling_pathways_from_reactome(
    reactome_data: dict,
    uniprot_id: str,
    gene_name: str,
) -> list[SignalingPathway]:
    """Convert Reactome fetch results into SignalingPathway models."""
    pathways = []
    seen: set[str] = set()

    for raw in reactome_data.get("signaling_pathways", []):
        pid = raw.get("stId") or raw.get("dbId") or ""
        if not pid or pid in seen:
            continue
        seen.add(str(pid))

        name = raw.get("displayName") or raw.get("name") or str(pid)
        category = classify_pathway(name)
        if category == "disease":
            continue  # handled separately

        pathways.append(SignalingPathway(
            pathway_id=str(pid),
            pathway_name=name,
            hierarchy=[],
            protein_role=gene_name,
            activates=[],
            inhibits=[],
            activated_by=[],
            inhibited_by=[],
            diagram_url=f"https://reactome.org/PathwayBrowser/#/{pid}",
            diagram_image_path=None,
        ))
    return pathways


def build_disease_pathways_from_reactome(
    reactome_data: dict,
    kegg_data: dict,
    gene_name: str,
) -> list[DiseasePathway]:
    """Convert Reactome + KEGG disease data into DiseasePathway models."""
    pathways = []
    seen: set[str] = set()

    # From Reactome disease_pathways
    for raw in reactome_data.get("disease_pathways", []):
        pid = str(raw.get("stId") or raw.get("dbId") or "")
        if not pid or pid in seen:
            continue
        seen.add(pid)
        name = raw.get("displayName") or raw.get("name") or pid
        pathways.append(DiseasePathway(
            disease_name=name,
            omim_id=None,
            pathway_id=pid,
            pathway_source="Reactome",
            cascade=[],
            diagram_url=f"https://reactome.org/PathwayBrowser/#/{pid}",
            diagram_image_path=None,
        ))

    # From KEGG disease IDs
    for did in kegg_data.get("disease_pathway_ids", []):
        clean = did.replace("ds:", "").replace("path:", "")
        if clean in seen:
            continue
        seen.add(clean)
        pathways.append(DiseasePathway(
            disease_name=clean,
            omim_id=None,
            pathway_id=clean,
            pathway_source="KEGG",
            cascade=[],
            diagram_url=f"https://www.genome.jp/dbget-bin/www_bget?ds:{clean}",
            diagram_image_path=None,
        ))

    return pathways


# ---------------------------------------------------------------------------
# PharmGKB pathway assembly
# ---------------------------------------------------------------------------

def build_pharmacogenomic_pathways_from_pharmgkb(
    pharmgkb_data: dict,
) -> list[PharmacogenomicPathway]:
    """Convert PharmGKB pathway data into PharmacogenomicPathway models."""
    pathways = []
    for raw in pharmgkb_data.get("pathways", []):
        pid = raw.get("id") or ""
        name = raw.get("name") or pid
        drugs = [
            c.get("name") or ""
            for c in raw.get("relatedChemicals", []) or []
        ]
        genes = [
            g.get("symbol") or g.get("name") or ""
            for g in raw.get("relatedGenes", []) or []
        ]
        pathways.append(PharmacogenomicPathway(
            pathway_id=pid,
            pathway_name=name,
            drugs_involved=[d for d in drugs if d],
            genes_involved=[g for g in genes if g],
            diagram_url=f"https://www.pharmgkb.org/pathway/{pid}",
            diagram_image_path=None,
        ))
    return pathways


# ---------------------------------------------------------------------------
# STRING interaction assembly
# ---------------------------------------------------------------------------

def build_protein_interactions_from_string(
    string_data: dict,
    min_score: float = 0.4,
) -> list[ProteinInteraction]:
    """Convert STRING interaction data into ProteinInteraction models."""
    interactions = []
    for item in string_data.get("interactions", []):
        score = item.get("score", 0)
        if isinstance(score, (int, float)):
            score_f = float(score) / 1000.0
        else:
            score_f = 0.0
        if score_f < min_score:
            continue
        partner = item.get("preferredName_B") or item.get("stringId_B") or ""
        if not partner:
            continue
        exp_score = item.get("escore", 0) or 0
        interactions.append(ProteinInteraction(
            partner_protein=partner,
            partner_gene=partner,
            interaction_type="functional",
            confidence_score=round(score_f, 3),
            experimentally_confirmed=float(exp_score) > 0,
        ))
    return interactions
