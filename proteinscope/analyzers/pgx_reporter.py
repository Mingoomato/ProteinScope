"""PGx Reporter — merge PharmGKB + DGIdb + ChEMBL data into structured models.

Synthesizes drug interaction data from multiple sources into DrugInteraction
and PGxVariant Pydantic models, applying optional evidence level filtering.
"""

from __future__ import annotations

from typing import Optional

from core.models import DrugInteraction, PGxVariant

# PharmGKB evidence level ordering (1A is highest confidence)
_LEVEL_ORDER = {"1A": 1, "1B": 2, "2A": 3, "2B": 4, "3": 5, "4": 6}


def _level_passes(level: str, min_level: str) -> bool:
    """Return True if the evidence level meets or exceeds min_level."""
    level_norm = level.upper().strip()
    min_norm = min_level.upper().strip()
    return _LEVEL_ORDER.get(level_norm, 99) <= _LEVEL_ORDER.get(min_norm, 99)


def build_drug_interactions(
    pharmgkb_data: dict,
    dgidb_interactions: list[dict],
    chembl_data: dict,
    min_pgx_evidence: str = "4",
) -> list[DrugInteraction]:
    """Build DrugInteraction models from all three drug-data sources.

    Args:
        pharmgkb_data:    Output of fetch_all_pharmgkb().
        dgidb_interactions: Output of fetch_dgidb_interactions().
        chembl_data:      Output of fetch_all_chembl().
        min_pgx_evidence: Minimum PharmGKB evidence level to include (default '4' = all).

    Returns:
        Deduplicated list of DrugInteraction models.
    """
    drug_map: dict[str, DrugInteraction] = {}

    # ── PharmGKB drug interactions ──────────────────────────────────────────
    for item in pharmgkb_data.get("drug_interactions", []):
        name = (item.get("name") or item.get("symbol") or "").strip()
        if not name:
            continue
        key = name.lower()
        drug_map[key] = DrugInteraction(
            drug_name=name,
            drug_id=item.get("id", ""),
            interaction_type=item.get("type", "unknown"),
            mechanism=None,
            clinical_significance=None,
            fda_label_mentioned=False,
            sources=["PharmGKB"],
        )

    # ── DGIdb interactions ───────────────────────────────────────────────────
    for item in dgidb_interactions:
        drug = item.get("drug") or {}
        name = (drug.get("name") or "").strip()
        if not name:
            continue
        key = name.lower()
        types = item.get("interactionTypes") or []
        itype = types[0].get("type", "unknown") if types else "unknown"
        if key in drug_map:
            existing = drug_map[key]
            if "DGIdb" not in existing.sources:
                existing.sources.append("DGIdb")
            if not existing.interaction_type or existing.interaction_type == "unknown":
                existing.interaction_type = itype
        else:
            drug_map[key] = DrugInteraction(
                drug_name=name,
                drug_id=drug.get("conceptId", ""),
                interaction_type=itype,
                mechanism=None,
                clinical_significance=None,
                fda_label_mentioned=False,
                sources=["DGIdb"],
            )

    # ── ChEMBL mechanisms ────────────────────────────────────────────────────
    for mech in chembl_data.get("mechanisms", []):
        name = (mech.get("molecule_name") or mech.get("molecule_chembl_id") or "").strip()
        if not name:
            continue
        key = name.lower()
        mechanism_text = mech.get("mechanism_of_action") or mech.get("action_type") or ""
        if key in drug_map:
            existing = drug_map[key]
            if "ChEMBL" not in existing.sources:
                existing.sources.append("ChEMBL")
            if mechanism_text and not existing.mechanism:
                existing.mechanism = mechanism_text
        else:
            drug_map[key] = DrugInteraction(
                drug_name=name,
                drug_id=mech.get("molecule_chembl_id", ""),
                interaction_type=mech.get("action_type", "unknown"),
                mechanism=mechanism_text or None,
                clinical_significance=None,
                fda_label_mentioned=False,
                sources=["ChEMBL"],
            )

    # ── Mark FDA label mentions ──────────────────────────────────────────────
    fda_drug_names = set()
    for label in pharmgkb_data.get("fda_labels", []):
        for chem in label.get("relatedChemicals", []):
            n = (chem.get("name") or "").strip().lower()
            if n:
                fda_drug_names.add(n)

    for key, di in drug_map.items():
        if key in fda_drug_names:
            di.fda_label_mentioned = True

    return list(drug_map.values())


def build_pgx_variants(
    pharmgkb_data: dict,
    gene_symbol: str,
    min_pgx_evidence: str = "4",
) -> list[PGxVariant]:
    """Build PGxVariant models from PharmGKB clinical annotation data.

    Args:
        pharmgkb_data:    Output of fetch_all_pharmgkb().
        gene_symbol:      HGNC gene symbol.
        min_pgx_evidence: Minimum PharmGKB evidence level to include.

    Returns:
        Filtered list of PGxVariant models.
    """
    variants: list[PGxVariant] = []
    for annotation in pharmgkb_data.get("pgx_variants", []):
        level = str(annotation.get("evidenceLevel") or annotation.get("level") or "4")
        if not _level_passes(level, min_pgx_evidence):
            continue

        # Extract variant rsID
        variant_id = ""
        for loc in annotation.get("location", {}).get("variants", []) or []:
            variant_id = loc.get("name") or loc.get("id") or ""
            break

        # Extract drug name
        drug_name = ""
        for chem in annotation.get("relatedChemicals", []) or []:
            drug_name = chem.get("name") or ""
            break

        # Extract phenotype
        phenotype = annotation.get("phenotype") or annotation.get("sentence") or ""
        if isinstance(phenotype, dict):
            phenotype = phenotype.get("term") or ""

        # Extract PubMed IDs
        pmids = [
            str(pub.get("pmid") or pub.get("id") or "")
            for pub in annotation.get("evidence", []) or []
            if pub.get("pmid") or pub.get("id")
        ]

        variants.append(PGxVariant(
            variant_id=variant_id or "unknown",
            gene_symbol=gene_symbol,
            drug_name=drug_name,
            phenotype=phenotype or "not specified",
            evidence_level=level,
            population=annotation.get("population"),
            pubmed_ids=pmids,
        ))
    return variants
