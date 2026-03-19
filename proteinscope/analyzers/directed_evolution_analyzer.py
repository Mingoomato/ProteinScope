"""Directed Evolution Design Analyzer.

Given a protein sequence and variant fitness scores (e.g. from ESM-2 ΔLL or
DMS assays), identifies hotspot residues and designs a suite of directed
evolution strategies with library size and round recommendations.

Background:
  Directed evolution is an iterative laboratory process alternating
  random/semi-random mutagenesis with functional selection to improve or
  re-programme enzyme/protein activity.  Rational selection of hotspot
  residues and appropriate mutagenesis strategy dramatically increases the
  probability of success.

  Key methods:
    Error-prone PCR — introduces random point mutations across the full gene;
      library sizes ~10^6; suitable when no structural rationale exists.
    Site-saturation mutagenesis (SSM) — exhaustively samples all 20 aa at
      chosen positions; library = 20^N for N positions; best for <5 hotspots.
    DNA shuffling (Stemmer) — recombines homologous sequences from multiple
      parents; ideal for multi-property optimisation of larger proteins.

Key citations:
  # Hotspot identification for directed evolution: Reetz MT 2006 Angew Chem
  # doi:10.1002/anie.200503017
  # Error-prone PCR: Cirino PC 2003 Methods Enzymol
  # doi:10.1016/S0076-6879(03)88003-6
  # DNA shuffling: Stemmer WP 1994 Nature doi:10.1038/370389a0
  # ESM-2 ΔLL as fitness proxy: Meier J 2021 bioRxiv
  # doi:10.1101/2021.07.09.450648
"""

from __future__ import annotations

import json as _json
import logging
from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel, Field

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class EvolutionStrategy(BaseModel):
    """A single directed evolution strategy recommendation."""

    method: str                             # "error-prone PCR" | "DNA shuffling" | "site-saturation mutagenesis"
    target_residues: List[int]
    estimated_library_size: int
    rounds_recommended: int
    selection_criteria: str
    expected_fold_improvement: Optional[str] = None    # e.g. "2-5x"
    citation: str = ""


class DirectedEvolutionPlan(BaseModel):
    """Full directed evolution design plan for a gene/protein."""

    gene: str
    hotspot_residues: List[int]
    strategies: List[EvolutionStrategy]
    recommended_strategy: str
    recommended_rationale: str
    gemini_rationale: str = ""
    timestamp: datetime


# ---------------------------------------------------------------------------
# Hotspot identification
# ---------------------------------------------------------------------------

def _identify_hotspots(
    variant_fitness_scores: list,
    binding_sites: list,
    active_sites: list,
) -> list[int]:
    """Identify hotspot residues from variant fitness scores.

    A position is a candidate hotspot if its delta_log_likelihood (ΔLL) > -1.0
    (i.e., mutations are reasonably tolerated at that position) AND the position
    overlaps with a known binding or active site.  If no overlap is found with
    functional sites, the top 10 most tolerant positions are returned as hotspots.

    # Hotspot identification for directed evolution: Reetz MT 2006 Angew Chem
    # doi:10.1002/anie.200503017
    # ESM-2 ΔLL as fitness proxy: Meier J 2021 bioRxiv doi:10.1101/2021.07.09.450648

    Tolerant threshold: ΔLL > -1.0
    # Threshold from: Hsu C 2022 Nat Biotechnol doi:10.1038/s41587-021-01122-5

    Args:
        variant_fitness_scores: List of dicts or objects with "position" and
                                 "delta_log_likelihood" (or "score") keys.
        binding_sites:          List of ints or dicts with "position" key.
        active_sites:           List of ints or dicts with "position" key.

    Returns:
        Sorted list of hotspot residue positions (1-based integers).
    """
    def _pos(obj) -> Optional[int]:
        if isinstance(obj, int):
            return obj
        if isinstance(obj, dict):
            p = obj.get("position") or obj.get("pos")
            return int(p) if p is not None else None
        p = getattr(obj, "position", None) or getattr(obj, "pos", None)
        return int(p) if p is not None else None

    def _dll(obj) -> float:
        if isinstance(obj, dict):
            v = obj.get("delta_log_likelihood") or obj.get("dll") or obj.get("score")
            return float(v) if v is not None else float("-inf")
        for attr in ("delta_log_likelihood", "dll", "score"):
            val = getattr(obj, attr, None)
            if val is not None:
                return float(val)
        return float("-inf")

    # Build set of functional site positions
    functional_positions: set[int] = set()
    for s in (binding_sites or []) + (active_sites or []):
        p = _pos(s)
        if p is not None:
            functional_positions.add(p)

    # Filter to tolerant positions (ΔLL > -1.0)
    # Threshold: Hsu C 2022 Nat Biotechnol doi:10.1038/s41587-021-01122-5
    tolerant: list[tuple[int, float]] = []
    for vf in variant_fitness_scores or []:
        pos = _pos(vf)
        dll = _dll(vf)
        if pos is None:
            continue
        if dll > -1.0:
            tolerant.append((pos, dll))

    # Sort by ΔLL descending (most tolerant first)
    tolerant.sort(key=lambda t: t[1], reverse=True)

    # Prefer positions overlapping functional sites
    if functional_positions:
        overlap = [pos for pos, _ in tolerant if pos in functional_positions]
        if overlap:
            return sorted(overlap)

    # Fallback: top 10 most tolerant positions
    return sorted(pos for pos, _ in tolerant[:10])


# ---------------------------------------------------------------------------
# Strategy design
# ---------------------------------------------------------------------------

def _design_strategies(
    hotspot_residues: list[int],
    sequence: str,
    active_sites: list,
) -> list[EvolutionStrategy]:
    """Propose three directed evolution strategies for the identified hotspots.

    Always returns exactly three strategies:
      a) Error-prone PCR targeting all hotspots.
      b) Site-saturation mutagenesis at the top 3 hotspot positions.
      c) DNA shuffling (if sequence > 100 aa).

    # Error-prone PCR: Cirino PC 2003 Methods Enzymol
    # doi:10.1016/S0076-6879(03)88003-6
    # DNA shuffling: Stemmer WP 1994 Nature doi:10.1038/370389a0
    # Site-saturation mutagenesis: Reetz MT 2008 Angew Chem
    # doi:10.1002/anie.200705157

    Args:
        hotspot_residues: Sorted list of hotspot positions.
        sequence:         Single-letter amino acid sequence.
        active_sites:     Raw active-site list for selection criteria text.

    Returns:
        List of three EvolutionStrategy objects.
    """
    strategies: list[EvolutionStrategy] = []
    top3_hotspots = hotspot_residues[:3] if hotspot_residues else []

    # ── a) Error-prone PCR ───────────────────────────────────────────────────
    strategies.append(EvolutionStrategy(
        method="error-prone PCR",
        target_residues=hotspot_residues,
        # Library size ~10^6 is standard for epPCR with Taq + MnCl2
        # Citation: Cirino PC 2003 Methods Enzymol doi:10.1016/S0076-6879(03)88003-6
        estimated_library_size=1_000_000,
        rounds_recommended=3,
        selection_criteria="functional assay / colorimetric screen",
        expected_fold_improvement="2-10x per round",
        citation=(
            "Cirino PC 2003 Methods Enzymol doi:10.1016/S0076-6879(03)88003-6"
        ),
    ))

    # ── b) Site-saturation mutagenesis at top 3 hotspots ────────────────────
    # Library size = 20^N for N positions (NNK codon degeneracy reduces this in
    # practice; 20^3 = 8000 is the theoretical maximum for 3 positions)
    # Citation: Reetz MT 2008 Angew Chem doi:10.1002/anie.200705157
    ssm_library = 20 ** max(len(top3_hotspots), 1)
    strategies.append(EvolutionStrategy(
        method="site-saturation mutagenesis",
        target_residues=top3_hotspots,
        estimated_library_size=ssm_library,
        rounds_recommended=2,
        selection_criteria=(
            "high-throughput plate assay at target positions; "
            "use NNK degenerate codons to reduce stop codon frequency"
        ),
        expected_fold_improvement="5-50x (position-specific optimisation)",
        citation=(
            "Reetz MT 2008 Angew Chem doi:10.1002/anie.200705157"
        ),
    ))

    # ── c) DNA shuffling (for sequences > 100 aa) ────────────────────────────
    # Citation: Stemmer WP 1994 Nature doi:10.1038/370389a0
    if len(sequence) > 100:
        strategies.append(EvolutionStrategy(
            method="DNA shuffling",
            target_residues=hotspot_residues,
            # DNase I fragmentation + random reassembly; typical library ~10^4
            estimated_library_size=10_000,
            rounds_recommended=4,
            selection_criteria=(
                "multi-property functional screen; suitable for recombining "
                "beneficial mutations from multiple epPCR variants"
            ),
            expected_fold_improvement="10-100x (cumulative beneficial recombination)",
            citation=(
                "Stemmer WP 1994 Nature doi:10.1038/370389a0"
            ),
        ))
    else:
        # Short sequence: substitute second round of epPCR with increased mutation rate
        strategies.append(EvolutionStrategy(
            method="DNA shuffling",
            target_residues=hotspot_residues,
            estimated_library_size=10_000,
            rounds_recommended=4,
            selection_criteria=(
                "Note: sequence < 100 aa — consider synthetic gene shuffling "
                "or StEP (staggered extension process) as alternative to "
                "standard DNase I shuffling"
            ),
            expected_fold_improvement="5-30x",
            citation=(
                "Stemmer WP 1994 Nature doi:10.1038/370389a0; "
                "Zhao H 1998 Nat Biotechnol doi:10.1038/nbt0398-258"
            ),
        ))

    return strategies


# ---------------------------------------------------------------------------
# Recommended strategy selection
# ---------------------------------------------------------------------------

def _select_recommended_strategy(
    hotspot_residues: list[int],
    active_sites: list,
) -> tuple[str, str]:
    """Select the preferred strategy and provide a one-sentence rationale.

    Decision rule:
      - If hotspot count < 5 AND any hotspot overlaps with an active site →
        site-saturation mutagenesis (focused, exhaustive sampling).
      - Otherwise → error-prone PCR (broader coverage).

    Args:
        hotspot_residues: Identified hotspot positions.
        active_sites:     Raw active-site objects/dicts.

    Returns:
        Tuple of (strategy_name, rationale_sentence).
    """
    def _pos(obj) -> Optional[int]:
        if isinstance(obj, int):
            return obj
        if isinstance(obj, dict):
            p = obj.get("position") or obj.get("pos")
            return int(p) if p is not None else None
        p = getattr(obj, "position", None) or getattr(obj, "pos", None)
        return int(p) if p is not None else None

    active_positions = set()
    for s in (active_sites or []):
        p = _pos(s)
        if p is not None:
            active_positions.add(p)

    hotspot_set   = set(hotspot_residues)
    active_overlap = hotspot_set & active_positions

    if len(hotspot_residues) < 5 and active_overlap:
        return (
            "site-saturation mutagenesis",
            (
                f"Fewer than 5 hotspot residues identified ({len(hotspot_residues)}) "
                f"with direct active-site overlap at position(s) "
                f"{sorted(active_overlap)} — exhaustive SSM maximises coverage "
                "of beneficial substitutions within a manageable library size."
            ),
        )
    return (
        "error-prone PCR",
        (
            f"{len(hotspot_residues)} hotspot(s) identified "
            f"{'without active-site overlap' if not active_overlap else 'spanning multiple functional regions'} "
            "— error-prone PCR provides broad sequence space exploration "
            "with a large library suitable for colorimetric or growth selection."
        ),
    )


# ---------------------------------------------------------------------------
# Gemini synthesis
# ---------------------------------------------------------------------------

async def _synthesize_evolution(
    gene: str,
    hotspot_residues: list[int],
    strategies: list[EvolutionStrategy],
    recommended: str,
    rationale: str,
) -> str:
    """Ask Gemini for detailed experimental design guidance.

    Args:
        gene:             Gene symbol.
        hotspot_residues: Identified hotspot positions.
        strategies:       Proposed evolution strategies.
        recommended:      Name of the recommended strategy.
        rationale:        Rule-based rationale sentence.

    Returns:
        Gemini narrative string; empty string on failure.
    """
    try:
        from core.gemini_interpreter import _call

        strat_json = _json.dumps(
            [
                {
                    "method":                   s.method,
                    "target_residues":          s.target_residues,
                    "estimated_library_size":   s.estimated_library_size,
                    "rounds_recommended":       s.rounds_recommended,
                    "selection_criteria":       s.selection_criteria,
                    "expected_fold_improvement": s.expected_fold_improvement,
                }
                for s in strategies
            ],
            indent=2,
        )

        prompt = (
            f"You are an expert in directed protein evolution. "
            f"For gene {gene}, the following hotspot residues have been identified: "
            f"{hotspot_residues}.\n\n"
            f"Proposed strategies:\n{strat_json}\n\n"
            f"Recommended strategy: {recommended}\n"
            f"Rationale: {rationale}\n\n"
            "Provide a detailed experimental design narrative covering:\n"
            "1. Specific mutagenesis protocol parameters for the recommended strategy "
            "(mutation rate, primer design for SSM, or DNase I conditions for shuffling).\n"
            "2. Recommended selection/screening assay type and throughput requirements.\n"
            "3. Expected improvement trajectory over the proposed rounds.\n"
            "4. Key quality-control checkpoints between rounds.\n"
            "5. Convergence criteria — when to stop iterating.\n\n"
            "Write as a concise expert laboratory protocol narrative (2-4 paragraphs). "
            "No markdown headers, no bullet points."
        )

        raw = await _call(prompt)
        return raw.strip() if raw else ""
    except Exception as exc:
        _log.debug("Gemini directed-evolution synthesis failed for %s: %s", gene, exc)
        return ""


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_directed_evolution_design(
    gene: str,
    sequence: str,
    variant_fitness_scores: list,
    binding_sites: list,
    active_sites: list,
    step_cb=None,
) -> Optional[DirectedEvolutionPlan]:
    """Design a directed evolution plan for a protein.

    Args:
        gene:                   Gene symbol, e.g. "TEM1".
        sequence:               Single-letter amino acid sequence.
        variant_fitness_scores: List of dicts/objects with "position" and
                                 "delta_log_likelihood" (or "score") fields.
        binding_sites:          List of binding-site positions (int or dict
                                 with "position" key).
        active_sites:           List of active-site positions (int or dict
                                 with "position" key).
        step_cb:                Optional async progress callback (receives str).

    Returns:
        DirectedEvolutionPlan, or None if a fatal error occurs.
    """
    async def _step(msg: str) -> None:
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    try:
        # ── 1. Identify hotspot residues ─────────────────────────────────────
        await _step("[1/4] Identifying hotspot residues...")
        hotspot_residues = _identify_hotspots(
            variant_fitness_scores,
            binding_sites,
            active_sites,
        )
        _log.debug(
            "Directed evolution hotspots for %s: %d positions",
            gene, len(hotspot_residues),
        )

        # ── 2. Design evolution strategies ───────────────────────────────────
        await _step("[2/4] Designing evolution strategies...")
        strategies = _design_strategies(hotspot_residues, sequence, active_sites)

        # ── 3. Select recommended strategy ───────────────────────────────────
        await _step("[3/4] Selecting recommended strategy...")
        recommended_strategy, recommended_rationale = _select_recommended_strategy(
            hotspot_residues, active_sites
        )

        # ── 4. Gemini synthesis ───────────────────────────────────────────────
        await _step("[4/4] Gemini synthesis...")
        gemini_rationale = await _synthesize_evolution(
            gene,
            hotspot_residues,
            strategies,
            recommended_strategy,
            recommended_rationale,
        )

        return DirectedEvolutionPlan(
            gene=gene,
            hotspot_residues=hotspot_residues,
            strategies=strategies,
            recommended_strategy=recommended_strategy,
            recommended_rationale=recommended_rationale,
            gemini_rationale=gemini_rationale,
            timestamp=datetime.utcnow(),
        )

    except Exception as exc:
        _log.debug("run_directed_evolution_design failed for %s: %s", gene, exc)
        return None
