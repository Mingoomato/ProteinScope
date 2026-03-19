"""IDP / LLPS Detection Analyzer.

Background — Intrinsically Disordered Proteins (IDPs):
─────────────────────────────────────────────────────────────────────────────
Intrinsically disordered proteins (IDPs) and intrinsically disordered regions
(IDRs) lack a stable three-dimensional structure under physiological conditions.
They exist as a conformational ensemble and are defined by:

  Per-residue disorder score (IUPred3):
    - Score > 0.5 for ≥20 consecutive residues → IDR
    - Long disordered regions (>30 residues) often serve as interaction hubs,
      linkers, or undergo disorder-to-order transitions upon binding

  Liquid–Liquid Phase Separation (LLPS) / FuzDrop:
    - IDRs enriched in aromatic (Tyr, Phe) or charged residues can undergo LLPS
    - FuzDrop predicts per-residue LLPS propensity; score > 0.6 → droplet-forming
    - LLPS underlies formation of membrane-less organelles (stress granules,
      P-bodies, nucleolus, transcriptional condensates)

  AlphaFold pLDDT integration:
    - pLDDT < 50: highly disordered, no stable fold
    - pLDDT 50–70: low-confidence, likely flexible linker or IDR
    - pLDDT > 70: predicted folded domain

Therapeutic relevance:
  - PROTACs / molecular glues: degrade disordered oncoproteins (e.g. c-MYC, AR-V7)
  - Stapled peptides: constrain disordered binding interfaces
  - Phase-separation disruptors: inhibit LLPS condensate formation (e.g. 1,6-hexanediol,
    small-molecule condensate inhibitors targeting FUS, TDP-43, hnRNPA1)
  - Molecular chaperones: prevent aberrant LLPS in neurodegeneration (TDP-43, FUS, tau)

This analyzer:
  1. Parses IUPred3 per-residue disorder scores
  2. Identifies IDRs (score > 0.5 for ≥20 residues)
  3. Parses FuzDrop LLPS propensity scores
  4. Correlates with AlphaFold pLDDT confidence
  5. Runs Gemini IDP functional synthesis
  6. Assembles IDP analysis report
"""

from __future__ import annotations

import asyncio
import json as _json
import logging
from typing import List, Optional

from pydantic import BaseModel, Field

from core.evidence import EvidenceGrade, DataProvenance

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

class IDPRegion(BaseModel):
    start: int
    end: int
    length: int
    region_type: str  # "disordered" | "llps_driver"
    mean_disorder_score: float
    mean_llps_score: float = 0.0
    provenance: Optional[DataProvenance] = None


class IDPAnalysis(BaseModel):
    gene: str
    sequence_length: int
    fraction_disordered: float
    fraction_llps_prone: float
    idr_regions: List[IDPRegion] = Field(default_factory=list)
    llps_overall_score: float = 0.0
    therapeutic_strategies: List[str] = Field(default_factory=list)
    functional_implications: str = ""
    gemini_synthesis: str = ""
    provenance: Optional[DataProvenance] = None


# ---------------------------------------------------------------------------
# Region detection helpers
# ---------------------------------------------------------------------------

def _find_idr_regions(scores: list[float], threshold: float = 0.5, min_len: int = 20) -> list[dict]:
    """Identify IDR regions from per-residue disorder scores."""
    regions = []
    in_region = False
    start = 0
    for i, sc in enumerate(scores):
        if sc > threshold and not in_region:
            in_region = True
            start = i
        elif sc <= threshold and in_region:
            if (i - start) >= min_len:
                regions.append({"start": start, "end": i - 1})
            in_region = False
    if in_region and (len(scores) - start) >= min_len:
        regions.append({"start": start, "end": len(scores) - 1})
    return regions


def _mean_slice(scores: list[float], start: int, end: int) -> float:
    """Mean of scores[start:end+1], returns 0.0 if empty."""
    if not scores or end < start:
        return 0.0
    sl = scores[start : end + 1]
    return round(sum(sl) / len(sl), 4) if sl else 0.0


# ---------------------------------------------------------------------------
# Gemini synthesis
# ---------------------------------------------------------------------------

async def _synthesize_idp(
    gene: str,
    fraction_disordered: float,
    idr_regions: list[IDPRegion],
    llps_score: float,
    af_plddt: Optional[float],
) -> tuple[str, list[str], str]:
    """Call Gemini for IDP functional synthesis.

    Returns:
        (functional_implications, therapeutic_strategies, gemini_synthesis)
    """
    try:
        from core.gemini_interpreter import _call

        idr_list = "; ".join(
            f"aa {r.start}-{r.end} (len={r.length}, type={r.region_type})"
            for r in idr_regions[:8]
        ) or "none detected"

        prompt = (
            f"Gene: {gene}\n"
            f"Fraction disordered: {fraction_disordered:.1%}\n"
            f"IDR regions: {idr_list}\n"
            f"LLPS score: {llps_score:.2f}\n"
            f"AlphaFold pLDDT: {af_plddt if af_plddt is not None else 'N/A'}\n\n"
            "You are an expert in intrinsically disordered proteins. "
            "Based on this IDP data:\n"
            "1. What are the likely functional roles of the disordered regions?\n"
            "2. What IDP-specific therapeutic strategies are relevant "
            "(PROTACs, stapled peptides, phase-separation disruptors, "
            "small-molecule chaperones)?\n"
            "3. What are key experimental approaches to validate these predictions?\n\n"
            "Return ONLY raw JSON (no markdown):\n"
            '{"functional_implications": "...", '
            '"therapeutic_strategies": ["...", "..."], '
            '"synthesis": "2-3 sentence summary"}'
        )

        raw = await _call(prompt)
        if not raw:
            return "AI synthesis unavailable.", [], "AI synthesis unavailable."

        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)

        functional_implications = str(data.get("functional_implications", "")).strip()
        strategies_raw = data.get("therapeutic_strategies", [])
        therapeutic_strategies = (
            [str(x).strip() for x in strategies_raw if x]
            if isinstance(strategies_raw, list)
            else []
        )
        gemini_synthesis = str(data.get("synthesis", "")).strip()

        return functional_implications, therapeutic_strategies, gemini_synthesis

    except Exception as exc:
        _log.debug("IDP Gemini synthesis failed: %s", exc)
        return "AI synthesis failed.", [], "AI synthesis failed."


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_idp_analysis(
    gene: str,
    sequence: str,
    iupred_data: dict,
    fuzdrop_data: dict,
    af_plddt: Optional[float] = None,
    step_cb=None,
) -> IDPAnalysis:
    """Run IDP / LLPS analysis pipeline.

    Args:
        gene:         Gene symbol (e.g. "FUS", "TDP43").
        sequence:     Canonical amino acid sequence.
        iupred_data:  Output of fetch_iupred_scores() — may be empty dict.
        fuzdrop_data: Output of fetch_fuzdrop() — may be empty dict.
        af_plddt:     AlphaFold overall pLDDT score (0-100), or None.
        step_cb:      Optional async progress callback(msg: str).

    Returns:
        IDPAnalysis — always returns; never raises.
    """

    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    seq_len = len(sequence) if sequence else 0

    # Graceful early exit if no data at all
    if not iupred_data and not fuzdrop_data:
        return IDPAnalysis(
            gene=gene.upper(),
            sequence_length=seq_len,
            fraction_disordered=0.0,
            fraction_llps_prone=0.0,
        )

    # ── 1. Parse IUPred3 disorder scores ────────────────────────────────────
    await _step("[1/6] Parsing IUPred3 disorder scores…")
    disorder_scores: list[float] = []
    if isinstance(iupred_data, dict):
        raw_scores = iupred_data.get("scores", [])
        disorder_scores = [float(s) for s in raw_scores if s is not None]

    # ── 2. Identify IDRs ─────────────────────────────────────────────────────
    await _step("[2/6] Identifying intrinsically disordered regions (IDRs)…")
    idr_raw = _find_idr_regions(disorder_scores, threshold=0.5, min_len=20)

    idr_provenance = DataProvenance(
        source="IUPred3 (iupred3.elte.hu)",
        evidence_grade=EvidenceGrade.COMPUTATIONAL,
        method="Long disorder prediction (score > 0.5 for ≥20 residues)",
        scientific_caveat=(
            "IUPred3 per-residue score; experimental IDR boundaries may differ"
        ),
    )

    idr_regions: list[IDPRegion] = []
    for reg in idr_raw:
        s, e = reg["start"], reg["end"]
        idr_regions.append(IDPRegion(
            start=s + 1,       # 1-based for output
            end=e + 1,
            length=e - s + 1,
            region_type="disordered",
            mean_disorder_score=_mean_slice(disorder_scores, s, e),
            provenance=idr_provenance,
        ))

    fraction_disordered = (
        sum(r.length for r in idr_regions) / seq_len if seq_len > 0 else 0.0
    )

    # ── 3. Parse FuzDrop LLPS scores ─────────────────────────────────────────
    await _step("[3/6] Parsing FuzDrop LLPS propensity scores…")
    per_residue_llps: list[float] = []
    llps_overall: float = 0.0
    droplet_regions: list[dict] = []

    if isinstance(fuzdrop_data, dict):
        per_residue_llps = [float(s) for s in fuzdrop_data.get("per_residue_scores", []) if s is not None]
        llps_overall = float(fuzdrop_data.get("llps_score", 0.0))
        droplet_regions = fuzdrop_data.get("droplet_regions", [])

    # Annotate LLPS scores onto existing IDR regions
    fuzdrop_provenance = DataProvenance(
        source="FuzDrop (fuzdrop.bio.unipd.it)",
        evidence_grade=EvidenceGrade.COMPUTATIONAL,
        method="Per-residue LLPS propensity; droplet-promoting score > 0.6",
        scientific_caveat=(
            "FuzDrop is a sequence-based predictor; in-vitro LLPS assay required for validation"
        ),
    )

    for region in idr_regions:
        s0 = region.start - 1   # back to 0-based
        e0 = region.end - 1
        region.mean_llps_score = _mean_slice(per_residue_llps, s0, e0)

    # Add LLPS-specific driver regions (droplet_regions not already in IDR list)
    existing_starts = {r.start for r in idr_regions}
    for dr in droplet_regions:
        dr_start = int(dr.get("start", 0))
        dr_end = int(dr.get("end", 0))
        if dr_start not in existing_starts and dr_end > dr_start:
            s0 = dr_start - 1
            e0 = dr_end - 1
            idr_regions.append(IDPRegion(
                start=dr_start,
                end=dr_end,
                length=dr_end - dr_start + 1,
                region_type="llps_driver",
                mean_disorder_score=_mean_slice(disorder_scores, s0, e0),
                mean_llps_score=_mean_slice(per_residue_llps, s0, e0),
                provenance=fuzdrop_provenance,
            ))

    # Fraction LLPS-prone residues (per_residue score > 0.6)
    llps_prone_count = sum(1 for sc in per_residue_llps if sc > 0.6)
    fraction_llps_prone = llps_prone_count / seq_len if seq_len > 0 else 0.0

    # ── 4. Correlate with AlphaFold pLDDT ───────────────────────────────────
    await _step("[4/6] Correlating with AlphaFold pLDDT confidence…")
    # If pLDDT is very low and no IUPred scores, treat entire sequence as disordered
    if af_plddt is not None and af_plddt < 50 and not idr_regions and seq_len > 0:
        idr_regions.append(IDPRegion(
            start=1,
            end=seq_len,
            length=seq_len,
            region_type="disordered",
            mean_disorder_score=1.0,
            mean_llps_score=llps_overall,
            provenance=DataProvenance(
                source="AlphaFold EBI v4 (inferred)",
                evidence_grade=EvidenceGrade.COMPUTATIONAL,
                confidence_interval=f"pLDDT={af_plddt:.1f}",
                scientific_caveat="Disorder inferred from very low pLDDT (<50); IUPred3 data unavailable",
            ),
        ))
        fraction_disordered = 1.0

    # ── 5. Gemini IDP functional synthesis ───────────────────────────────────
    await _step("[5/6] Running Gemini IDP functional synthesis…")
    functional_implications, therapeutic_strategies, gemini_synthesis = (
        await _synthesize_idp(gene, fraction_disordered, idr_regions, llps_overall, af_plddt)
    )

    # ── 6. Assemble report ───────────────────────────────────────────────────
    await _step("[6/6] Assembling IDP analysis report…")
    overall_provenance = DataProvenance(
        source="IUPred3 + FuzDrop + AlphaFold + Gemini 2.5 Pro",
        evidence_grade=EvidenceGrade.COMPUTATIONAL,
        scientific_caveat=(
            "Disorder and LLPS predictions are sequence-based; "
            "experimental validation (NMR, single-molecule FRET, in-vitro LLPS assay) required"
        ),
    )

    return IDPAnalysis(
        gene=gene.upper(),
        sequence_length=seq_len,
        fraction_disordered=round(fraction_disordered, 4),
        fraction_llps_prone=round(fraction_llps_prone, 4),
        idr_regions=idr_regions,
        llps_overall_score=round(llps_overall, 4),
        therapeutic_strategies=therapeutic_strategies,
        functional_implications=functional_implications,
        gemini_synthesis=gemini_synthesis,
        provenance=overall_provenance,
    )
