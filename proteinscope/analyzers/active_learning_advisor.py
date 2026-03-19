"""Active learning advisor for protein engineering mutation campaigns.

Implements a Gaussian Process-inspired Upper Confidence Bound (UCB) acquisition
function to recommend which mutations to measure next in a directed evolution or
DMS campaign.

Gaussian Process active learning for protein fitness:
Frazier PI 2018 arXiv 1807.02811

Usage::

    from analyzers.active_learning_advisor import run_active_learning_advice

    recommendation = await run_active_learning_advice(
        gene="BRCA1",
        sequence="MDLSALRVEEVQ...",
        variant_fitness_scores=scored_variants,
        mavedb_data=await fetch_mavedb_scores("BRCA1"),
    )
    for mut in recommendation.recommended_next:
        print(mut.position, mut.alt_aa, mut.acquisition_score)
"""

from __future__ import annotations

import math
from datetime import datetime
from typing import Dict, List, Optional

from pydantic import BaseModel, Field

from core.evidence import DataProvenance, EvidenceGrade


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


class RecommendedMutation(BaseModel):
    position: int
    ref_aa: str
    alt_aa: str
    expected_improvement: float
    uncertainty: float
    acquisition_score: float    # UCB = expected_improvement + 1.5 * uncertainty
    rationale: str
    provenance: Optional[DataProvenance] = None


class MutationRecommendation(BaseModel):
    gene: str
    recommended_next: List[RecommendedMutation] = Field(default_factory=list)   # top 5
    mavedb_overlap_count: int
    gemini_strategy: str = ""
    timestamp: datetime


# ---------------------------------------------------------------------------
# Amino acids considered as alternatives for UCB recommendation
# ---------------------------------------------------------------------------
# Canonical substitution candidates (chemically diverse, commonly mutagenized)
_ALT_AAS = list("AVLIGSTC")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _get_attr(obj, key: str, default=None):
    """Safely retrieve attribute or dict key."""
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def _parse_position(variant_str: str) -> Optional[int]:
    """Extract 1-based position integer from a variant string like 'A123V'."""
    import re
    m = re.search(r"(\d+)", str(variant_str))
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return None
    return None


def _parse_ref_aa(variant_str: str) -> str:
    """Extract reference amino acid from a variant string like 'A123V'."""
    import re
    m = re.match(r"^([A-Za-z])", str(variant_str))
    if m:
        return m.group(1).upper()
    return "X"


def _parse_alt_aa(variant_str: str) -> str:
    """Extract alt amino acid from a variant string like 'A123V'."""
    import re
    m = re.search(r"[A-Za-z]$", str(variant_str))
    if m:
        return m.group(0).upper()
    return "X"


def _build_position_landscape(
    variant_fitness_scores: list,
) -> Dict[int, Dict[str, list]]:
    """Build {position: {alt_aa: [delta_ll, ...]}} from fitness scores.

    Also stores ref_aa per position.
    """
    landscape: Dict[int, Dict[str, list]] = {}
    for item in variant_fitness_scores:
        try:
            variant_str = _get_attr(item, "variant", None)
            dll = _get_attr(item, "delta_log_likelihood", None)
            if variant_str is None or dll is None:
                continue
            pos = _parse_position(variant_str)
            if pos is None:
                continue
            alt_aa = _parse_alt_aa(variant_str)
            ref_aa = _parse_ref_aa(variant_str)
            dll = float(dll)

            if pos not in landscape:
                landscape[pos] = {"_ref": ref_aa, "_dlls": [dll]}
            else:
                landscape[pos]["_dlls"].append(dll)

            if alt_aa not in landscape[pos]:
                landscape[pos][alt_aa] = []
            landscape[pos][alt_aa].append(dll)
        except Exception:
            continue
    return landscape


def _mean(values: list[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def _std(values: list[float]) -> float:
    if len(values) < 2:
        return 1.5  # high uncertainty for singletons
    m = _mean(values)
    variance = sum((v - m) ** 2 for v in values) / len(values)
    return math.sqrt(variance) if variance > 0 else 0.1


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


async def run_active_learning_advice(
    gene: str,
    sequence: str,
    variant_fitness_scores: list,
    mavedb_data: list,
    step_cb=None,
) -> Optional[MutationRecommendation]:
    """Recommend next mutations to measure using a UCB acquisition function.

    Implements a GP-surrogate UCB strategy:
      UCB(i, a) = mean_ΔLL(i, a) + 1.5 * σ(i)

    This balances exploitation (high mean fitness) with exploration (high
    uncertainty regions not yet sampled).

    Citation: Gaussian Process active learning for protein fitness:
    Frazier PI 2018 arXiv 1807.02811

    Args:
        gene:                   Gene symbol (e.g. "BRCA1").
        sequence:               Full amino acid sequence string.
        variant_fitness_scores: List of VariantScore-like objects or dicts with
                                variant + delta_log_likelihood keys.
        mavedb_data:            Output from fetch_mavedb_scores().
        step_cb:                Optional async progress callback (receives str).

    Returns:
        MutationRecommendation on success, None on unrecoverable failure.
    """

    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    try:
        seq_upper = sequence.upper() if sequence else ""

        # ── Step 1: Build observed mutation landscape ────────────────────────
        await _step("[1/4] Building observed mutation landscape...")

        landscape = _build_position_landscape(variant_fitness_scores)

        # ── Step 2: Estimate GP uncertainty per position ─────────────────────
        await _step("[2/4] Estimating GP uncertainty...")

        # Citation: Gaussian Process active learning for protein fitness:
        # Frazier PI 2018 arXiv 1807.02811
        # For positions with measurements: uncertainty = std(ΔLL values)
        #   - singleton → 1.5 (high uncertainty)
        # For unmeasured positions: uncertainty = 2.0 (uninformative prior)

        position_uncertainty: Dict[int, float] = {}
        for pos, data in landscape.items():
            dlls = data.get("_dlls", [])
            if len(dlls) == 0:
                position_uncertainty[pos] = 2.0
            elif len(dlls) == 1:
                position_uncertainty[pos] = 1.5  # singleton — high uncertainty
            else:
                position_uncertainty[pos] = max(_std(dlls), 0.1)

        # Also consider all positions in sequence as candidates
        # (positions not yet measured get uncertainty = 2.0)
        all_positions = set(range(1, len(seq_upper) + 1))
        measured_positions = set(landscape.keys())
        unmeasured_positions = all_positions - measured_positions

        # ── Step 3: Compute UCB acquisition function ─────────────────────────
        await _step("[3/4] Computing acquisition function (UCB)...")

        # Citation: UCB acquisition function:
        # Frazier PI 2018 arXiv 1807.02811
        # UCB(i, alt_aa) = mean_ΔLL(i, alt_aa) + 1.5 * σ(i)
        UCB_BETA = 1.5  # exploration-exploitation tradeoff

        acquisition_candidates: list[tuple[int, str, str, float, float, float]] = []
        # (position, ref_aa, alt_aa, mean_dll, uncertainty, ucb_score)

        # Candidates from measured positions with high uncertainty or alt_aa options
        for pos, data in landscape.items():
            ref_aa = data.get("_ref", "X")
            if isinstance(ref_aa, list):
                ref_aa = ref_aa[0] if ref_aa else "X"
            uncertainty = position_uncertainty.get(pos, 1.5)

            for alt_aa in _ALT_AAS:
                if alt_aa == ref_aa:
                    continue
                alt_dlls = data.get(alt_aa, [])
                if alt_dlls:
                    mean_dll = _mean(alt_dlls)
                else:
                    # Not yet measured for this alt_aa at this position
                    mean_dll = 0.0
                ucb = mean_dll + UCB_BETA * uncertainty
                acquisition_candidates.append(
                    (pos, ref_aa, alt_aa, mean_dll, uncertainty, ucb)
                )

        # Candidates from unmeasured positions (top 20% of sequence by position index)
        # Limit to a reasonable number of unmeasured positions to avoid combinatorial explosion
        unmeasured_sample = sorted(unmeasured_positions)[:50]
        for pos in unmeasured_sample:
            seq_idx = pos - 1
            ref_aa = seq_upper[seq_idx] if 0 <= seq_idx < len(seq_upper) else "X"
            uncertainty = 2.0  # uninformative prior
            for alt_aa in _ALT_AAS:
                if alt_aa == ref_aa:
                    continue
                mean_dll = 0.0  # prior mean — no data
                ucb = mean_dll + UCB_BETA * uncertainty
                acquisition_candidates.append(
                    (pos, ref_aa, alt_aa, mean_dll, uncertainty, ucb)
                )

        # Sort by UCB descending
        acquisition_candidates.sort(key=lambda x: x[5], reverse=True)

        # Deduplicate by position — take highest UCB per position, then pick top 5
        seen_positions: set[int] = set()
        top_candidates: list[tuple[int, str, str, float, float, float]] = []
        for cand in acquisition_candidates:
            pos = cand[0]
            if pos not in seen_positions:
                seen_positions.add(pos)
                top_candidates.append(cand)
            if len(top_candidates) >= 5:
                break

        prov = DataProvenance(
            source="GP-UCB active learning surrogate (ΔLL landscape)",
            evidence_grade=EvidenceGrade.COMPUTATIONAL,
            scientific_caveat=(
                "UCB scores are derived from a GP approximation over ESM-2 ΔLL values; "
                "actual fitness may differ — experimental validation required."
            ),
            method="UCB = mean_ΔLL + 1.5 * σ (Frazier 2018 arXiv:1807.02811)",
        )

        recommended: list[RecommendedMutation] = []
        for pos, ref_aa, alt_aa, mean_dll, uncertainty, ucb in top_candidates:
            measured_tag = "measured" if pos in measured_positions else "not yet measured"
            rationale = (
                f"Position {pos} ({ref_aa}→{alt_aa}): UCB={ucb:.3f}, "
                f"mean ΔLL={mean_dll:.3f}, uncertainty σ={uncertainty:.3f}. "
                f"Position is {measured_tag}. "
                f"High UCB indicates {'exploitation of beneficial signal' if mean_dll > 0 else 'exploration of uncertain region'}."
            )
            recommended.append(
                RecommendedMutation(
                    position=pos,
                    ref_aa=ref_aa,
                    alt_aa=alt_aa,
                    expected_improvement=round(mean_dll, 4),
                    uncertainty=round(uncertainty, 4),
                    acquisition_score=round(ucb, 4),
                    rationale=rationale,
                    provenance=prov,
                )
            )

        # ── Step 4: Cross-reference MaveDB ───────────────────────────────────
        await _step("[4/4] Cross-referencing MaveDB...")

        mavedb_overlap_count = 0
        for item in mavedb_data:
            try:
                urn = _get_attr(item, "urn", None)
                target_gene = _get_attr(item, "targetGene", None)
                if urn or target_gene:
                    mavedb_overlap_count += 1
            except Exception:
                continue

        # Gemini synthesis
        gemini_text = ""
        try:
            from core.gemini_interpreter import _call

            top_lines = "\n".join(
                f"  {m.position} ({m.ref_aa}→{m.alt_aa}): UCB={m.acquisition_score:.3f}, "
                f"mean ΔLL={m.expected_improvement:.3f}, σ={m.uncertainty:.3f}"
                for m in recommended
            ) if recommended else "  No candidates computed."

            prompt = (
                f"Gene: {gene}\n"
                f"MaveDB score sets available: {mavedb_overlap_count}\n"
                f"Top 5 recommended mutations (UCB active learning):\n{top_lines}\n\n"
                "You are a protein engineering expert specializing in directed evolution "
                "and active learning for fitness landscape exploration.\n"
                "1. Provide a prioritized experimental strategy for these mutations.\n"
                "2. Which mutations should be tested first and why (epistasis risk, "
                "structural context, therapeutic relevance)?\n"
                "3. What assay format (yeast display, PACE, cell viability) would be "
                "most efficient for this gene?\n\n"
                "Return a concise 3-4 sentence experimental prioritization strategy."
            )

            gemini_text = await _call(prompt) or ""
        except Exception:
            gemini_text = ""

        return MutationRecommendation(
            gene=gene,
            recommended_next=recommended,
            mavedb_overlap_count=mavedb_overlap_count,
            gemini_strategy=gemini_text.strip(),
            timestamp=datetime.utcnow(),
        )

    except Exception:
        return None
