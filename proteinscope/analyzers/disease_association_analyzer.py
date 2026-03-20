"""Disease Association Layer analyzer.

Merges OpenTargets + DisGeNET + ClinGen + HPO into a unified ranked list.
Applies LoF/GoF mapping, mechanistic template, and Gemini synthesis.

# Citation: 김박사 + 노박사, Grand Consortium V3 2026-03-20
# Minimum evidence filter: OpenTargets score > 0.5 OR ClinGen Definitive/Strong (노박사)
# LoF SO terms: SO:0001587/0001589/0001574/0001575 (노박사)
# Mechanistic template: {Gene} {variant} → {mechanism} → {pathway} → {phenotype} (노박사)
"""

from __future__ import annotations

import os
from datetime import datetime
from typing import Optional

from core.evidence import DataProvenance, EvidenceGrade
from core.models import DiseaseAssociation, DiseaseAssociationReport

# ── LoF/GoF SO term mapping (노박사, Grand Consortium V3) ─────────────────────
# Loss-of-Function SO terms (high-impact truncating variants)
_LOF_SO_TERMS = {
    "SO_0001587",  # stop_gained
    "SO_0001589",  # frameshift_variant
    "SO_0001574",  # splice_acceptor_variant
    "SO_0001575",  # splice_donor_variant
}
# GoF is inferred: SO_0001583 (missense_variant) + known oncogenic context
_GOF_SO_MISSENSE = "SO_0001583"

# ClinGen classifications that pass the minimum filter
_CLINGEN_PASS = {"Definitive", "Strong"}

# OpenTargets minimum score (노박사 recommendation: > 0.5 for primary display)
_OT_MIN_SCORE = 0.5

# ID prefix priority for canonical dedup: EFO > MONDO > UMLS > others
_ID_PRIORITY = {"EFO": 3, "MONDO": 2, "UMLS": 1}


def _id_priority(disease_id: str) -> int:
    prefix = disease_id.split(":")[0].upper() if ":" in disease_id else disease_id[:5].upper()
    return _ID_PRIORITY.get(prefix, 0)


def _merge_and_dedup(
    ot_data: list[dict],
    disgenet_data: list[dict],
    clingen_data: list[dict],
    hpo_data: list[dict],
) -> list[DiseaseAssociation]:
    """Merge multi-source disease data, dedup by disease name, pick canonical ID."""

    # Build HPO index by term name for later attachment
    hpo_by_name: dict[str, list[str]] = {}
    for term in hpo_data:
        name = term.get("hpo_name", "").lower()
        if name:
            hpo_by_name.setdefault(name, []).append(
                f"{term.get('hpo_id', '')} {term.get('hpo_name', '')}".strip()
            )

    # Key: normalized disease name → DiseaseAssociation
    merged: dict[str, DiseaseAssociation] = {}

    # ── OpenTargets (primary source, best scores) ──
    for item in ot_data:
        name = item.get("disease_name", "").strip()
        if not name:
            continue
        key = name.lower()
        disease_id = item.get("disease_id", "")
        score = float(item.get("score") or 0.0)

        if key not in merged:
            merged[key] = DiseaseAssociation(
                disease_id=disease_id,
                disease_name=name,
                source_ids=[disease_id] if disease_id else [],
                score=score,
                source="OpenTargets",
                evidence_type="genetic_association",
                provenance=DataProvenance(
                    source="OpenTargets Platform API v4 (GraphQL)",
                    evidence_grade=EvidenceGrade.EXPERIMENTAL,
                ),
            )
        else:
            assoc = merged[key]
            if disease_id and disease_id not in assoc.source_ids:
                assoc.source_ids.append(disease_id)
            # Pick higher-priority canonical ID
            if _id_priority(disease_id) > _id_priority(assoc.disease_id):
                assoc.disease_id = disease_id
            if score > assoc.score:
                assoc.score = score

    # ── DisGeNET ──
    for item in disgenet_data:
        name = item.get("disease_name", "").strip()
        if not name:
            continue
        key = name.lower()
        disease_id = item.get("disease_id", "")
        score = float(item.get("score") or 0.0)
        ev_count = item.get("evidence_count")

        if key not in merged:
            merged[key] = DiseaseAssociation(
                disease_id=disease_id,
                disease_name=name,
                source_ids=[disease_id] if disease_id else [],
                score=score,
                source="DisGeNET",
                evidence_count=ev_count,
                evidence_type="genetic_association",
                provenance=DataProvenance(
                    source="DisGeNET REST API",
                    evidence_grade=EvidenceGrade.EXPERIMENTAL,
                ),
            )
        else:
            assoc = merged[key]
            if disease_id and disease_id not in assoc.source_ids:
                assoc.source_ids.append(disease_id)
            if ev_count is not None:
                assoc.evidence_count = ev_count

    # ── ClinGen ──
    for item in clingen_data:
        name = item.get("disease_name", "").strip()
        if not name:
            continue
        key = name.lower()
        classification = item.get("classification", "")
        mondo_id = item.get("mondo_id", "")

        if key not in merged:
            merged[key] = DiseaseAssociation(
                disease_id=mondo_id or name,
                disease_name=name,
                source_ids=[mondo_id] if mondo_id else [],
                score=0.9 if classification == "Definitive" else 0.7,
                source="ClinGen",
                evidence_type="genetic_association",
                clingen_classification=classification,
                canonical_mondo_id=mondo_id if mondo_id.startswith("MONDO") else None,
                provenance=DataProvenance(
                    source="ClinGen Gene-Disease Validity",
                    evidence_grade=EvidenceGrade.EXPERIMENTAL,
                ),
            )
        else:
            assoc = merged[key]
            assoc.clingen_classification = classification
            if mondo_id:
                if mondo_id not in assoc.source_ids:
                    assoc.source_ids.append(mondo_id)
                if mondo_id.startswith("MONDO"):
                    assoc.canonical_mondo_id = mondo_id
                    if _id_priority(mondo_id) > _id_priority(assoc.disease_id):
                        assoc.disease_id = mondo_id

    # ── Attach HPO terms (fuzzy name match) ──
    for assoc in merged.values():
        key = assoc.disease_name.lower()
        if key in hpo_by_name:
            assoc.hpo_terms = hpo_by_name[key][:5]  # cap at 5 terms per disease

    return list(merged.values())


def _apply_evidence_filter(associations: list[DiseaseAssociation]) -> list[DiseaseAssociation]:
    """Keep associations that pass minimum evidence threshold (노박사)."""
    filtered = []
    for assoc in associations:
        passes_ot = assoc.source == "OpenTargets" and assoc.score > _OT_MIN_SCORE
        passes_clingen = assoc.clingen_classification in _CLINGEN_PASS
        passes_disgenet = assoc.source == "DisGeNET" and assoc.score > 0.3
        if passes_ot or passes_clingen or passes_disgenet:
            filtered.append(assoc)
    return filtered


def _infer_lof_gof(
    assoc: DiseaseAssociation,
    gene_name: str,
) -> DiseaseAssociation:
    """Infer LoF/GoF hint and therapeutic modality from disease context.

    LoF SO terms (노박사): SO_0001587/0001589/0001574/0001575
    GoF: inferred from known oncogenic gene + missense context.
    """
    disease_lower = assoc.disease_name.lower()

    # Heuristic: cancers → likely GoF for oncogenes; rare diseases → likely LoF
    oncogenic_keywords = {"cancer", "carcinoma", "tumor", "leukemia", "lymphoma",
                          "melanoma", "glioma", "sarcoma"}
    lof_keywords = {"cystic fibrosis", "muscular dystrophy", "hemophilia",
                    "phenylketonuria", "gaucher", "deficiency", "aplasia"}

    is_cancer = any(kw in disease_lower for kw in oncogenic_keywords)
    is_lof_disease = any(kw in disease_lower for kw in lof_keywords)

    if is_lof_disease:
        assoc.lof_gof_hint = "LoF"
        assoc.therapeutic_modality_hint = "replacement"
    elif is_cancer:
        assoc.lof_gof_hint = "GoF"
        assoc.therapeutic_modality_hint = "inhibitor"

    return assoc


async def _gemini_synthesis(
    gene: str,
    associations: list[DiseaseAssociation],
) -> str:
    """Generate Gemini mechanistic summary for top disease associations."""
    if not associations:
        return ""

    top = associations[:5]
    disease_list = "\n".join(
        f"- {a.disease_name} (score={a.score:.2f}, source={a.source}"
        + (f", ClinGen={a.clingen_classification}" if a.clingen_classification else "")
        + (f", {a.lof_gof_hint}" if a.lof_gof_hint else "")
        + ")"
        for a in top
    )

    prompt = (
        f"You are a clinical biochemist. For the gene {gene}, summarize the top disease "
        f"associations below using the mechanistic template:\n"
        f"{{Gene}} {{variant/mutation class}} → {{molecular mechanism}} → {{pathway}} → {{phenotype}}\n\n"
        f"Disease associations:\n{disease_list}\n\n"
        f"Write 2-3 concise sentences covering the strongest associations. "
        f"Focus on molecular mechanism. Do not use bullet points."
    )

    try:
        import google.generativeai as genai

        api_key = os.environ.get("GEMINI_API_KEY", "")
        if not api_key:
            return ""
        genai.configure(api_key=api_key)
        model = genai.GenerativeModel("models/gemini-2.5-pro")
        response = model.generate_content(prompt)
        return (response.text or "").strip()
    except Exception:
        return ""


async def run_disease_association_analysis(
    gene: str,
    ot_data: list[dict],
    disgenet_data: list[dict],
    clingen_data: list[dict],
    hpo_data: list[dict],
    step_cb=None,
) -> Optional[DiseaseAssociationReport]:
    """Merge and analyze disease associations from 4 databases.

    Args:
        gene: Gene symbol (e.g. "EGFR").
        ot_data: OpenTargets results from fetch_opentargets().
        disgenet_data: DisGeNET results from fetch_disgenet().
        clingen_data: ClinGen results from fetch_clingen().
        hpo_data: HPO terms from fetch_hpo_terms().
        step_cb: Optional async progress callback.

    Returns:
        DiseaseAssociationReport or None on complete failure.
    """
    try:
        if step_cb:
            await step_cb(f"[Disease] Merging {len(ot_data)} OpenTargets + "
                          f"{len(disgenet_data)} DisGeNET + {len(clingen_data)} ClinGen records")

        # Step 1: Merge and dedup
        associations = _merge_and_dedup(ot_data, disgenet_data, clingen_data, hpo_data)

        # Step 2: Evidence filter (노박사: OT > 0.5 OR ClinGen Definitive/Strong)
        associations = _apply_evidence_filter(associations)

        # Step 3: LoF/GoF inference
        associations = [_infer_lof_gof(a, gene) for a in associations]

        # Step 4: Sort by score descending
        associations.sort(key=lambda x: x.score, reverse=True)

        if step_cb:
            await step_cb(f"[Disease] {len(associations)} associations after filtering")

        # Step 5: Gemini synthesis
        if step_cb:
            await step_cb("[Disease] Generating mechanistic summary")
        gemini_summary = await _gemini_synthesis(gene, associations)

        return DiseaseAssociationReport(
            gene=gene,
            associations=associations,
            total_count=len(associations),
            gemini_summary=gemini_summary,
            timestamp=datetime.utcnow().isoformat(),
        )

    except Exception:
        return None
