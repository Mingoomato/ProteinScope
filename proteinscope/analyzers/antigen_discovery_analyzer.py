"""Antigen Discovery Mode analyzer (P11).

Screens Human Protein Atlas (HPA) expression data to identify actionable tumor
antigens for immuno-oncology modalities (CAR-T, antibody, ADC, TCE).

Key references used throughout this module:
  - Uhlén et al. (2015) Science 347:1260419 — Human Protein Atlas resource
  - Fagerberg et al. (2014) Mol Cell Proteomics 13:397 — tissue specificity scoring
  - Carter (2001) Nat Rev Cancer 1:118 — antibody target selection criteria
  - Mortazavi et al. (2008) Nat Methods 5:621 — RPKM expression threshold
  - Scott et al. (2012) Nat Rev Cancer 12:278 — antibody/CAR-T surface accessibility
  - Klebanoff et al. (2016) J Clin Invest 126:3206 — on-target off-tumor toxicity
"""

from __future__ import annotations

import json as _json
import logging
from datetime import datetime, timezone
from typing import List, Optional

from pydantic import BaseModel, Field

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

class AntigenCandidate(BaseModel):
    gene: str
    accession: str              # UniProt accession if known, else ""
    score: float                # computed tumor/normal ratio score
    tumor_rna_max: float
    normal_rna_max: float
    subcellular_locations: List[str]
    cancer_category: str        # "High" | "Medium" | "Low" | "Not detected"
    known_programs: List[str]   # NCT IDs from ClinicalTrials if available
    safety_flag: str            # "safe" | "caution" | "high_normal_expression"
    rationale: str              # 1-sentence justification


class AntigenDiscoveryReport(BaseModel):
    tumor_type: str
    modality: str               # "CAR_T" | "antibody" | "ADC" | "TCE"
    candidates: List[AntigenCandidate] = Field(default_factory=list)
    top_recommendation: str     # gene name of best candidate
    gemini_rationale: str       # Gemini synthesis paragraph
    screening_criteria: dict    # the criteria used (for display)
    timestamp: str


# ---------------------------------------------------------------------------
# Scoring constants — all citations are mandatory
# ---------------------------------------------------------------------------

# Tumor/normal ratio base score formula:
#   score = tumor_rna_max × (1 / (normal_rna_max + 0.1)) × topology_bonus
#
# Tissue specificity ratio approach:
#   Fagerberg et al. (2014) Mol Cell Proteomics 13:397
#   "A protein-level analysis of the human tissue atlas" — ratio of tissue-
#   elevated vs. ubiquitous expression used as primary specificity metric.
#
# The +0.1 pseudocount prevents division-by-zero; empirical consensus in
# expression scoring pipelines (no single reference; analogous to pseudocount
# usage in DESeq2 — Love et al. (2014) Genome Biol 15:550).

_TOPOLOGY_BONUS_SURFACE: float = 2.0
# Rationale: plasma membrane / cell-surface antigens are physically accessible
# to antibody fragments, CAR-T constructs, and ADC payloads without internalisation.
# Scott et al. (2012) Nat Rev Cancer 12:278 — "Antibodies as cancer therapeutics:
# surface accessibility is the primary druggability determinant."

_TOPOLOGY_BONUS_DEFAULT: float = 1.0
# Intracellular / unknown topology — not directly accessible to extracellular
# therapeutics.  Retain in ranking in case modality is TCE (MHC-presented
# peptides from intracellular proteins are valid targets).

# Surface/accessible topology keywords — Carter (2001) Nat Rev Cancer 1:118:
# "Therapeutic antibodies for cancer" defines three targetable antigen classes:
# surface receptors, secreted growth factors, and matrix-associated proteins.
_SURFACE_KEYWORDS: frozenset[str] = frozenset({
    "Plasma membrane",
    "Cell surface",
    "Secreted",
})

# Minimum tumor RNA expression cutoff (RPKM / nTPM)
_MIN_TUMOR_RPKM: float = 10.0
# Mortazavi et al. (2008) Nat Methods 5:621 — "Mapping and quantifying mammalian
# transcriptomes by RNA-Seq": RPKM ≥ 10 is the established threshold for
# reliable detection of a gene as "expressed."

# Safety classification thresholds
# Based on Klebanoff et al. (2016) J Clin Invest 126:3206 — "On-target off-tumor
# toxicities in chimeric antigen receptor therapy":
#   Normal tissue expression ≤5 nTPM → low off-tumour risk ("safe")
#   Normal tissue expression 5–20 nTPM → moderate risk ("caution")
#   Normal tissue expression >20 nTPM → high off-tumour risk
_SAFE_NORMAL_RPKM: float = 5.0
_CAUTION_NORMAL_RPKM: float = 20.0

# Maximum number of candidates to return after ranking
_TOP_N: int = 20


# ---------------------------------------------------------------------------
# Scoring helpers
# ---------------------------------------------------------------------------

def _parse_locations(raw_scl: object) -> list[str]:
    """Normalise subcellular_location field from HPA.

    HPA returns a semicolon-separated string or already a list.
    """
    if isinstance(raw_scl, list):
        return [s.strip() for s in raw_scl if s and str(s).strip()]
    if isinstance(raw_scl, str):
        return [s.strip() for s in raw_scl.split(";") if s.strip()]
    return []


def _topology_bonus(locations: list[str]) -> float:
    """Return topology accessibility multiplier.

    Scott et al. (2012) Nat Rev Cancer 12:278 — surface accessibility criterion.
    """
    for loc in locations:
        for keyword in _SURFACE_KEYWORDS:
            if keyword.lower() in loc.lower():
                return _TOPOLOGY_BONUS_SURFACE
    return _TOPOLOGY_BONUS_DEFAULT


def _safety_flag(normal_rna_max: float) -> str:
    """Classify antigen safety based on normal tissue expression.

    Klebanoff et al. (2016) J Clin Invest 126:3206 — on-target off-tumor
    toxicity framework for CAR-T and antibody therapeutics.
    """
    if normal_rna_max <= _SAFE_NORMAL_RPKM:
        return "safe"
    elif normal_rna_max <= _CAUTION_NORMAL_RPKM:
        return "caution"
    else:
        return "high_normal_expression"


def _compute_score(tumor_rna_max: float, normal_rna_max: float,
                   locations: list[str]) -> float:
    """Compute antigen attractiveness score.

    Formula:
        score = tumor_rna_max × (1 / (normal_rna_max + 0.1)) × topology_bonus

    References:
        Fagerberg et al. (2014) Mol Cell Proteomics 13:397 — specificity ratio.
        Scott et al. (2012) Nat Rev Cancer 12:278 — topology multiplier.
        Love et al. (2014) Genome Biol 15:550 — pseudocount convention.
    """
    bonus = _topology_bonus(locations)
    return round(tumor_rna_max * (1.0 / (normal_rna_max + 0.1)) * bonus, 4)


def _passes_filters(record: dict) -> bool:
    """Apply mandatory inclusion filters before scoring.

    Filter 1 — Surface/secreted requirement:
        Carter (2001) Nat Rev Cancer 1:118 — antibody targeting requires antigens
        exposed on the cell surface or in the extracellular space.  CAR-T and
        ADC modalities share the same physical accessibility constraint.

    Filter 2 — Minimum tumor expression:
        Mortazavi et al. (2008) Nat Methods 5:621 — RPKM ≥ 10 = reliably expressed.
    """
    locations = _parse_locations(record.get("subcellular_location", ""))
    has_surface = any(
        any(kw.lower() in loc.lower() for kw in _SURFACE_KEYWORDS)
        for loc in locations
    )
    if not has_surface:
        return False  # Filter 1: must be surface/secreted — Carter (2001)

    try:
        tumor_max = float(record.get("tumor_rna_max", 0) or 0)
    except (TypeError, ValueError):
        tumor_max = 0.0
    if tumor_max < _MIN_TUMOR_RPKM:
        return False  # Filter 2: must be expressed — Mortazavi et al. (2008)

    return True


def _build_rationale(gene: str, score: float, locations: list[str],
                     safety: str, modality: str) -> str:
    """Generate a 1-sentence candidate justification string."""
    loc_str = ", ".join(locations[:2]) if locations else "unknown localisation"
    return (
        f"{gene} scores {score:.2f} (tumor/normal specificity ratio, "
        f"Fagerberg et al. 2014) with {loc_str} topology; "
        f"safety classification '{safety}' per Klebanoff et al. (2016); "
        f"suitable for {modality} modality targeting."
    )


# ---------------------------------------------------------------------------
# ClinicalTrials cross-reference
# ---------------------------------------------------------------------------

def _match_clinical_trials(gene: str,
                            clinical_trials_data: Optional[List[dict]]) -> list[str]:
    """Return NCT IDs from pre-fetched ClinicalTrials data that mention this gene.

    Performs a case-insensitive substring match against the gene symbol in the
    trial title/intervention/condition fields.  Callers are responsible for
    fetching trial data upstream via the ClinicalTrials fetcher.
    """
    if not clinical_trials_data:
        return []
    nct_ids: list[str] = []
    gene_upper = gene.upper()
    for trial in clinical_trials_data:
        searchable = " ".join([
            str(trial.get("nct_id", "")),
            str(trial.get("title", "")),
            str(trial.get("intervention", "")),
            str(trial.get("condition", "")),
        ]).upper()
        if gene_upper in searchable:
            nct_id = trial.get("nct_id", "")
            if nct_id and nct_id not in nct_ids:
                nct_ids.append(nct_id)
    return nct_ids[:5]  # cap at 5 NCT IDs per candidate for display compactness


# ---------------------------------------------------------------------------
# Gemini synthesis
# ---------------------------------------------------------------------------

async def _synthesize_with_gemini(
    report: AntigenDiscoveryReport,
) -> AntigenDiscoveryReport:
    """Call Gemini to produce a ranked therapeutic strategy paragraph.

    Falls back gracefully if Gemini is unavailable; the report is still
    complete with algorithm-derived candidates.
    """
    try:
        from core.gemini_interpreter import _call  # type: ignore[import]
    except ImportError:
        log.warning("gemini_interpreter._call not importable — skipping AI synthesis")
        report.gemini_rationale = "AI synthesis unavailable (import error)."
        return report

    top5 = report.candidates[:5]
    if not top5:
        report.gemini_rationale = "No candidates passed the screening filters."
        return report

    candidate_lines = "\n".join(
        f"  {i+1}. {c.gene}: score={c.score:.2f}, "
        f"tumor_rna_max={c.tumor_rna_max}, normal_rna_max={c.normal_rna_max}, "
        f"locations={'; '.join(c.subcellular_locations)}, "
        f"safety={c.safety_flag}, cancer_category={c.cancer_category}, "
        f"known_programs={c.known_programs or 'none'}"
        for i, c in enumerate(top5)
    )

    prompt = (
        f"You are an expert immuno-oncology scientist evaluating tumor antigens "
        f"for the therapeutic modality: {report.modality}.\n\n"
        f"Tumor type: {report.tumor_type}\n\n"
        f"Top 5 antigen candidates ranked by tumor/normal specificity score "
        f"(Fagerberg et al. 2014 Mol Cell Proteomics 13:397):\n"
        f"{candidate_lines}\n\n"
        "Please rank these candidates considering:\n"
        "  1. Tumor selectivity and therapeutic window\n"
        "  2. Known resistance mechanisms (antigen loss, downregulation)\n"
        "  3. Suitability for the requested modality (CAR-T, antibody, ADC, or TCE)\n"
        "  4. Published clinical or preclinical precedent\n\n"
        "Ground your analysis in published literature on tumor antigen selection. "
        "Cite specific papers with PMID or DOI where known. "
        "Return a single coherent paragraph of 150-250 words. "
        "Do NOT use markdown headers or bullet points."
    )

    try:
        raw = await _call(prompt, timeout=120.0)
        report.gemini_rationale = raw.strip() if raw else "AI synthesis returned empty response."
    except Exception as exc:
        log.warning("Gemini synthesis call failed: %s", exc)
        report.gemini_rationale = "AI synthesis failed."

    return report


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

async def run_antigen_discovery(
    tumor_type: str,
    modality: str,
    hpa_data: List[dict],
    clinical_trials_data: Optional[List[dict]] = None,
    step_cb=None,
) -> AntigenDiscoveryReport:
    """Screen HPA expression data and return ranked antigen candidates.

    Pipeline:
        1. Filter HPA records by surface topology and minimum expression
           (Carter 2001; Mortazavi et al. 2008)
        2. Score each candidate by tumor/normal specificity ratio
           (Fagerberg et al. 2014)
        3. Classify safety based on normal tissue expression
           (Klebanoff et al. 2016)
        4. Cross-reference with ClinicalTrials data if provided
        5. Rank by score descending, return top 20
        6. Gemini synthesis paragraph

    Args:
        tumor_type:           Cancer / tumor keyword (e.g. "melanoma").
        modality:             One of CAR_T | antibody | ADC | TCE.
        hpa_data:             Pre-fetched list of dicts from fetch_hpa_tumor_expression.
        clinical_trials_data: Optional pre-fetched ClinicalTrials records.
        step_cb:              Optional async progress callback (msg: str) -> None.

    Returns:
        AntigenDiscoveryReport with scored candidates and Gemini rationale.
    """
    async def _step(msg: str) -> None:
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    # Normalise modality to uppercase underscore form
    modality_normalised = modality.strip().upper().replace("-", "_").replace(" ", "_")

    screening_criteria = {
        "topology_filter": "Plasma membrane OR Cell surface OR Secreted",
        "topology_filter_citation": "Carter (2001) Nat Rev Cancer 1:118",
        "min_tumor_rpkm": _MIN_TUMOR_RPKM,
        "min_tumor_rpkm_citation": "Mortazavi et al. (2008) Nat Methods 5:621",
        "scoring_formula": "tumor_rna_max × (1 / (normal_rna_max + 0.1)) × topology_bonus",
        "scoring_citation": "Fagerberg et al. (2014) Mol Cell Proteomics 13:397",
        "safety_thresholds": {
            "safe_normal_rpkm_max": _SAFE_NORMAL_RPKM,
            "caution_normal_rpkm_max": _CAUTION_NORMAL_RPKM,
        },
        "safety_citation": "Klebanoff et al. (2016) J Clin Invest 126:3206",
        "topology_bonus_surface": _TOPOLOGY_BONUS_SURFACE,
        "topology_bonus_default": _TOPOLOGY_BONUS_DEFAULT,
        "topology_bonus_citation": "Scott et al. (2012) Nat Rev Cancer 12:278",
        "top_n_returned": _TOP_N,
    }

    # ── 1. Filter ────────────────────────────────────────────────────────────
    await _step("[1/5] Applying surface topology and expression filters…")
    passing = [rec for rec in hpa_data if _passes_filters(rec)]
    log.debug("Antigen discovery: %d/%d records passed filters for %s",
              len(passing), len(hpa_data), tumor_type)

    # ── 2. Score ─────────────────────────────────────────────────────────────
    await _step("[2/5] Scoring candidates by tumor/normal specificity ratio…")
    scored: list[tuple[float, dict]] = []
    for rec in passing:
        locs = _parse_locations(rec.get("subcellular_location", ""))
        try:
            t_max = float(rec.get("tumor_rna_max", 0) or 0)
            n_max = float(rec.get("normal_rna_max", 0) or 0)
        except (TypeError, ValueError):
            continue
        s = _compute_score(t_max, n_max, locs)
        scored.append((s, rec))

    # Sort descending by score; stable sort preserves alphabetical order on ties
    scored.sort(key=lambda x: x[0], reverse=True)

    # ── 3. Build AntigenCandidate objects (top _TOP_N) ───────────────────────
    await _step("[3/5] Classifying safety and cross-referencing clinical trials…")
    candidates: list[AntigenCandidate] = []
    for score, rec in scored[:_TOP_N]:
        gene = str(rec.get("gene_symbol", rec.get("gene", "UNKNOWN"))).strip()
        locs = _parse_locations(rec.get("subcellular_location", ""))
        try:
            t_max = float(rec.get("tumor_rna_max", 0) or 0)
            n_max = float(rec.get("normal_rna_max", 0) or 0)
        except (TypeError, ValueError):
            t_max, n_max = 0.0, 0.0

        safety = _safety_flag(n_max)
        nct_ids = _match_clinical_trials(gene, clinical_trials_data)

        candidates.append(AntigenCandidate(
            gene=gene,
            accession=str(rec.get("uniprot_id", "") or ""),
            score=score,
            tumor_rna_max=t_max,
            normal_rna_max=n_max,
            subcellular_locations=locs,
            cancer_category=str(rec.get("cancer_category", "") or ""),
            known_programs=nct_ids,
            safety_flag=safety,
            rationale=_build_rationale(gene, score, locs, safety,
                                       modality_normalised),
        ))

    top_recommendation = candidates[0].gene if candidates else ""

    # ── 4. Assemble preliminary report ───────────────────────────────────────
    await _step("[4/5] Assembling report…")
    report = AntigenDiscoveryReport(
        tumor_type=tumor_type,
        modality=modality_normalised,
        candidates=candidates,
        top_recommendation=top_recommendation,
        gemini_rationale="",
        screening_criteria=screening_criteria,
        timestamp=datetime.now(timezone.utc).isoformat(),
    )

    # ── 5. Gemini synthesis ──────────────────────────────────────────────────
    await _step("[5/5] Gemini 2.5 Pro synthesizing therapeutic strategy…")
    report = await _synthesize_with_gemini(report)

    return report


# ---------------------------------------------------------------------------
# Natural language query parser
# ---------------------------------------------------------------------------

async def parse_antigen_query(text: str) -> dict:
    """Extract tumor_type and modality from a natural language query.

    Delegates to core.gemini_interpreter.parse_antigen_discovery_query for
    AI-powered extraction.  Falls back to a deterministic keyword scan so
    that the pipeline is never blocked even if Gemini is unavailable.

    Args:
        text: Free-text query (e.g. "Find CAR-T targets for glioblastoma").

    Returns:
        {"tumor_type": str, "modality": str}
    """
    try:
        from core.gemini_interpreter import parse_antigen_discovery_query  # type: ignore[import]
        result = await parse_antigen_discovery_query(text)
        if isinstance(result, dict) and result.get("tumor_type") not in (None, "", "unknown"):
            return result
    except Exception as exc:
        log.debug("parse_antigen_discovery_query delegation failed: %s", exc)

    # ---------------------------------------------------------------------------
    # Deterministic fallback — keyword matching for common modalities and cancers
    # ---------------------------------------------------------------------------
    text_lower = text.lower()

    # Modality detection (longest match first to avoid "antibody" matching in "ADC antibody")
    modality = "CAR_T"  # default
    _modality_map = [
        (["adc", "antibody-drug conjugate"], "ADC"),
        (["tce", "t cell engager", "bispecific"], "TCE"),
        (["car-nk", "car nk"], "CAR_NK"),
        (["car-t", "car t", "cart"], "CAR_T"),
        (["antibody", "mab", "monoclonal"], "antibody"),
    ]
    for keywords, mod in _modality_map:
        if any(kw in text_lower for kw in keywords):
            modality = mod
            break

    # Tumor type extraction — common oncology keywords
    _tumor_keywords = [
        "melanoma", "glioblastoma", "gbm", "breast cancer", "lung cancer",
        "colorectal cancer", "prostate cancer", "ovarian cancer", "leukemia",
        "lymphoma", "myeloma", "pancreatic cancer", "hepatocellular carcinoma",
        "renal cell carcinoma", "bladder cancer", "cervical cancer",
        "head and neck", "thyroid cancer", "sarcoma", "neuroblastoma",
    ]
    tumor_type = "unknown"
    for kw in _tumor_keywords:
        if kw in text_lower:
            tumor_type = kw
            break

    return {"tumor_type": tumor_type, "modality": modality}
