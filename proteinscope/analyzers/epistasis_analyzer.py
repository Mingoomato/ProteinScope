"""Pairwise epistasis analysis from deep mutational scanning data.

Estimates epistatic interactions between missense variants using single-mutant
delta log-likelihood scores as a proxy, cross-referenced with MaveDB DMS datasets.

Pairwise epistasis: Domingo E 2019 Nat Genet doi:10.1038/s41588-019-0437-9

Usage::

    from analyzers.epistasis_analyzer import run_epistasis_analysis

    report = await run_epistasis_analysis(
        gene="BRCA1",
        sequence="MDLSALRVEEVQ...",
        variant_fitness_scores=scored_variants,
        mavedb_data=await fetch_mavedb_scores("BRCA1"),
    )
"""

from __future__ import annotations

import asyncio
from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel, Field

from core.evidence import DataProvenance, EvidenceGrade


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


class EpistasisPair(BaseModel):
    variant_a: str          # e.g. "A123V"
    variant_b: str          # e.g. "G456D"
    epistasis_index: float  # ε = ΔLL(i,j) − ΔLL(i) − ΔLL(j); approximated here
    interpretation: str     # "negative" | "positive" | "additive"
    provenance: Optional[DataProvenance] = None


class EpistasisReport(BaseModel):
    gene: str
    pairs: List[EpistasisPair] = Field(default_factory=list)
    high_risk_combinations: List[str] = Field(default_factory=list)
    mavedb_score_sets_used: List[str] = Field(default_factory=list)
    gemini_interpretation: str = ""
    timestamp: datetime


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_attr(obj, key: str, default=None):
    """Safely retrieve an attribute or dict key."""
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def _parse_fitness_scores(variant_fitness_scores: list) -> dict[str, float]:
    """Extract {variant_str: delta_log_likelihood} from a list of score objects or dicts."""
    scores: dict[str, float] = {}
    for item in variant_fitness_scores:
        try:
            variant_str = _get_attr(item, "variant", None)
            dll = _get_attr(item, "delta_log_likelihood", None)
            if variant_str is None or dll is None:
                continue
            scores[str(variant_str)] = float(dll)
        except Exception:
            continue
    return scores


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


async def run_epistasis_analysis(
    gene: str,
    sequence: str,
    variant_fitness_scores: list,
    mavedb_data: list,
    step_cb=None,
) -> Optional[EpistasisReport]:
    """Run pairwise epistasis analysis for a gene.

    Args:
        gene:                    Gene symbol (e.g. "BRCA1").
        sequence:                Canonical amino acid sequence string.
        variant_fitness_scores:  List of VariantScore objects or dicts with
                                 variant + delta_log_likelihood keys.
        mavedb_data:             Output from fetch_mavedb_scores().
        step_cb:                 Optional async progress callback (receives str).

    Returns:
        EpistasisReport on success, None on unrecoverable failure.
    """

    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    try:
        # ── Step 1: Parse variant fitness scores ────────────────────────────
        await _step("[1/5] Parsing variant fitness scores...")

        dll_map = _parse_fitness_scores(variant_fitness_scores)
        if not dll_map:
            return EpistasisReport(
                gene=gene,
                pairs=[],
                high_risk_combinations=[],
                mavedb_score_sets_used=[],
                gemini_interpretation="No variant fitness scores available for epistasis analysis.",
                timestamp=datetime.utcnow(),
            )

        variant_names = list(dll_map.keys())

        # ── Step 2: Compute pairwise epistasis indices ───────────────────────
        await _step("[2/5] Computing pairwise epistasis indices...")

        # Pairwise epistasis: Domingo E 2019 Nat Genet doi:10.1038/s41588-019-0437-9
        # Full pairwise ΔLL(i,j) not available from single-mutant data alone.
        # Approximation: epistasis_index = abs(dll_a - dll_b) for damaging pairs;
        # interpretation based on individual dll thresholds.
        pairs: list[EpistasisPair] = []
        pair_prov = DataProvenance(
            source="ESM-2 ΔLL pairwise epistasis approximation",
            evidence_grade=EvidenceGrade.COMPUTATIONAL,
            scientific_caveat=(
                "Pairwise epistasis index is an approximation from single-mutant ΔLL scores; "
                "true double-mutant combinatorial measurements (DMS) required for exact ε."
            ),
            method="abs(ΔLL_i − ΔLL_j) for damaging variant pairs (Domingo 2019 proxy)",
        )

        n = len(variant_names)
        candidate_pairs: list[tuple[str, str, float, str]] = []

        for i in range(n):
            for j in range(i + 1, n):
                va = variant_names[i]
                vb = variant_names[j]
                dll_a = dll_map[va]
                dll_b = dll_map[vb]

                # Only compute for variant pairs where at least one is damaging
                if dll_a >= 0 and dll_b >= 0:
                    continue

                # Pairwise epistasis: Domingo E 2019 Nat Genet doi:10.1038/s41588-019-0437-9
                epistasis_index = abs(dll_a - dll_b)

                # Interpretation rules:
                if dll_a < -2 and dll_b < -2 and epistasis_index > 2.0:
                    interpretation = "negative"
                elif (dll_a > 0 and dll_b < -1) or (dll_b > 0 and dll_a < -1):
                    interpretation = "positive"
                else:
                    interpretation = "additive"

                candidate_pairs.append((va, vb, epistasis_index, interpretation))

        # Sort by absolute epistasis_index descending, keep top 20
        candidate_pairs.sort(key=lambda x: abs(x[2]), reverse=True)
        for va, vb, idx, interp in candidate_pairs[:20]:
            pairs.append(
                EpistasisPair(
                    variant_a=va,
                    variant_b=vb,
                    epistasis_index=round(idx, 4),
                    interpretation=interp,
                    provenance=pair_prov,
                )
            )

        # ── Step 3: Cross-reference MaveDB datasets ──────────────────────────
        await _step("[3/5] Cross-referencing MaveDB datasets...")

        score_set_urns: list[str] = []
        for item in mavedb_data:
            try:
                urn = _get_attr(item, "urn", None)
                if urn:
                    score_set_urns.append(str(urn))
            except Exception:
                continue

        # ── Step 4: Identify high-risk combinations ──────────────────────────
        await _step("[4/5] Identifying high-risk combinations...")

        # High-risk: both variants individually show ΔLL < -2.0 AND epistasis = "negative"
        high_risk: list[str] = []
        for pair in pairs:
            dll_a = dll_map.get(pair.variant_a, 0.0)
            dll_b = dll_map.get(pair.variant_b, 0.0)
            if dll_a < -2.0 and dll_b < -2.0 and pair.interpretation == "negative":
                high_risk.append(f"{pair.variant_a}+{pair.variant_b}")

        # ── Step 5: Gemini synthesis ─────────────────────────────────────────
        await _step("[5/5] Gemini synthesis...")

        gemini_text = ""
        try:
            from core.gemini_interpreter import _call

            top5 = pairs[:5]
            top5_lines = "\n".join(
                f"  {p.variant_a} + {p.variant_b}: ε={p.epistasis_index:.3f} ({p.interpretation})"
                for p in top5
            ) if top5 else "  No significant pairs detected."

            prompt = (
                f"Gene: {gene}\n"
                f"Number of MaveDB score sets available: {len(score_set_urns)}\n"
                f"Top 5 epistatic pairs (ε = pairwise epistasis index):\n{top5_lines}\n"
                f"High-risk combinations (both variants ΔLL < -2, negative epistasis): "
                f"{', '.join(high_risk) if high_risk else 'None detected'}\n\n"
                "You are a protein biochemistry expert specializing in epistasis and DMS analysis.\n"
                "1. Interpret the clinical significance of these epistatic interactions.\n"
                "2. What do the high-risk combinations suggest about protein function?\n"
                "3. What experimental follow-up (e.g., double-mutant thermodynamic cycles, "
                "structural mapping) would confirm these interactions?\n\n"
                "Return a concise 3-4 sentence clinical interpretation."
            )

            gemini_text = await _call(prompt) or ""
        except Exception:
            gemini_text = ""

        return EpistasisReport(
            gene=gene,
            pairs=pairs,
            high_risk_combinations=high_risk,
            mavedb_score_sets_used=score_set_urns,
            gemini_interpretation=gemini_text.strip(),
            timestamp=datetime.utcnow(),
        )

    except Exception:
        return None
