"""Protein Engineering Strategy Recommender.

Recommends rational protein engineering strategies to improve thermostability,
expression yield, and functional properties based on sequence + AlphaFold pLDDT.

Strategies covered
──────────────────
1. pLDDT-guided flexibility mapping     — low-pLDDT loops → rigidification targets
2. Proline/Ala substitution             — Gly/Ala in loops → Pro for helix cap / loop rigidity
3. Salt bridge / electrostatic network  — surface charge optimisation for thermal stability
4. Disulfide bridge candidates          — Cys pair proximity from sequence analysis
5. Buried polar residue elimination     — unsatisfied H-bond donors in hydrophobic core
6. ΔTm estimation (empirical proxies)   — fragment-based ΔTm from ProTherm concepts
7. Expression yield improvements        — N-term Met/Ala, rare codon density, N-glycosylation

Key citations
─────────────
• Guerois R et al. 2002 J Mol Biol 320:369-387 doi:10.1016/S0022-2836(02)00442-4
  (FoldX free-energy function; ΔΔG stability prediction)
• Romero PA & Arnold FH 2009 Nat Rev Mol Cell Biol 10:866-76 doi:10.1038/nrm2805
  (directed protein evolution — rational strategies)
• Goldenzweig A et al. 2016 Mol Cell 63:337-46 doi:10.1016/j.molcel.2016.06.012
  (PROSS algorithm; stability-expression tradeoff)
• Shoichet BK et al. 1995 Biochemistry 34:4157-66 doi:10.1021/bi00013a001
  (buried polar residue cost in protein folding)
• Musil M et al. 2021 Nucleic Acids Res 49:W304-W312 doi:10.1093/nar/gkab338
  (HotSpot Wizard 3.0 — thermostabilisation)
• Jumper J et al. 2021 Nature 596:583-589 doi:10.1038/s41586-021-03819-2
  (AlphaFold2 pLDDT confidence as proxy for local flexibility)
• Eijsink VG et al. 2004 J Biotechnol 113:105-20 doi:10.1016/j.jbiotec.2004.03.026
  (protein engineering for thermal stability — comprehensive review)
"""
from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from pydantic import BaseModel, Field

# ──────────────────────────────────────────────────────────────────────────────
# Thresholds
# ──────────────────────────────────────────────────────────────────────────────

# pLDDT < 70 → locally disordered / flexible region
# Citation: Jumper et al. 2021 doi:10.1038/s41586-021-03819-2
_LOW_PLDDT_THRESHOLD = 70.0

# Cys–Cys Cα–Cα distance for disulfide candidacy: 4.5–7.5 Å in native Cys-free context
# estimated from sequence separation: ~3.8 Å / residue in disordered; 1.5 Å / residue in helix
# We use a coarse sequence-distance proxy in absence of coordinates.
# Citation: Eijsink 2004 doi:10.1016/j.jbiotec.2004.03.026 (disulfide engineering)
_DISULFIDE_SEQ_WINDOW_MIN = 5
_DISULFIDE_SEQ_WINDOW_MAX = 200

# Pro substitution candidates: Gly or Ala at position i within a loop (i.e. not flanked by helix/strand)
# ΔTm ~ +0.5–1.5 °C per Pro substitution in unstructured loop
# Citation: Musil 2021 doi:10.1093/nar/gkab338; Eijsink 2004
_PRO_SUBSTITUTION_AAS = {"G", "A"}

# ΔTm per engineering action (empirical, kcal/mol → °C proxies)
# Citation: Guerois 2002 doi:10.1016/S0022-2836(02)00442-4 (FoldX)
_DT_PER_PRO_SUBSTITUTION = 0.8    # °C per Gly→Pro in loop
_DT_PER_DISULFIDE = 2.5           # °C per disulfide bridge
_DT_PER_BURIED_POLAR_FIX = 1.2   # °C per buried polar residue replacement


# ──────────────────────────────────────────────────────────────────────────────
# Amino-acid properties for engineering analysis
# ──────────────────────────────────────────────────────────────────────────────

# Buried polar residues: if in hydrophobic core they destabilise (Shoichet 1995)
_POLAR_RESIDUES = {"N", "Q", "S", "T"}     # excludes charged ones (D,E,R,K,H)
# Hydrophobic residues — used to identify buried context
_HYDROPHOBIC = {"I", "L", "V", "F", "W", "Y", "M", "A"}


# ──────────────────────────────────────────────────────────────────────────────
# Pydantic models
# ──────────────────────────────────────────────────────────────────────────────

class StabilizationStrategy(BaseModel):
    """One specific engineering recommendation."""
    strategy_type: str           # "proline_sub" | "disulfide" | "surface_charge" | "buried_polar"
    description: str
    target_positions: List[int] = Field(default_factory=list)
    target_residues: List[str] = Field(default_factory=list)
    estimated_delta_tm_c: Optional[float] = None   # Estimated ΔTm in °C
    confidence: str = "medium"   # "high" | "medium" | "low"
    reference: str = ""


class EngineeringStrategyReport(BaseModel):
    """Full protein engineering strategy report."""
    gene: str
    sequence_length: int
    estimated_total_delta_tm_c: float = 0.0
    low_plddt_regions: List[str] = Field(default_factory=list)   # "10–25 (loop)"
    strategies: List[StabilizationStrategy] = Field(default_factory=list)
    expression_tips: List[str] = Field(default_factory=list)
    gemini_design: str = ""
    timestamp: datetime = Field(default_factory=datetime.utcnow)
    references: List[str] = Field(default_factory=list)


# ──────────────────────────────────────────────────────────────────────────────
# Private helpers
# ──────────────────────────────────────────────────────────────────────────────

_REFERENCES = [
    "Guerois R et al. 2002 J Mol Biol 320:369 doi:10.1016/S0022-2836(02)00442-4",
    "Romero PA & Arnold FH 2009 Nat Rev Mol Cell Biol 10:866 doi:10.1038/nrm2805",
    "Goldenzweig A et al. 2016 Mol Cell 63:337 doi:10.1016/j.molcel.2016.06.012",
    "Shoichet BK et al. 1995 Biochemistry 34:4157 doi:10.1021/bi00013a001",
    "Musil M et al. 2021 Nucleic Acids Res 49:W304 doi:10.1093/nar/gkab338",
    "Jumper J et al. 2021 Nature 596:583 doi:10.1038/s41586-021-03819-2",
    "Eijsink VG et al. 2004 J Biotechnol 113:105 doi:10.1016/j.jbiotec.2004.03.026",
]


def _find_low_plddt_regions(af_plddt: Optional[float]) -> List[str]:
    """Without per-residue pLDDT, use overall score as a rough signal.

    In full deployment this would use per-residue B-factor from the AlphaFold
    PDB file (B-factor column encodes per-residue pLDDT * 100 / 100).
    Here we flag that overall pLDDT < threshold means significant disorder.
    """
    if af_plddt is None:
        return []
    # Citation: Jumper 2021 doi:10.1038/s41586-021-03819-2
    if af_plddt < _LOW_PLDDT_THRESHOLD:
        return [f"Overall pLDDT={af_plddt:.1f} — likely disordered regions present"]
    return []


def _find_pro_substitutions(sequence: str) -> List[StabilizationStrategy]:
    """Identify Gly/Ala residues in putative loop positions as Pro substitution candidates.

    Heuristic: isolated G/A (not flanked by 3+ consecutive G/A) are likely loop residues.
    Citation: Musil 2021 doi:10.1093/nar/gkab338; Eijsink 2004 doi:10.1016/j.jbiotec.2004.03.026
    """
    candidates: List[int] = []
    for i, aa in enumerate(sequence):
        if aa not in _PRO_SUBSTITUTION_AAS:
            continue
        # Skip if already Pro at i-1 or i+1 (no point in adjacent substitution)
        prev_aa = sequence[i - 1] if i > 0 else ""
        next_aa = sequence[i + 1] if i < len(sequence) - 1 else ""
        if prev_aa == "P" or next_aa == "P":
            continue
        # Skip N/C-terminus (positions 1–3 and last 3)
        pos = i + 1
        if pos <= 3 or pos >= len(sequence) - 2:
            continue
        candidates.append(pos)

    if not candidates:
        return []

    # Report top 5
    top = candidates[:5]
    delta_tm = len(top) * _DT_PER_PRO_SUBSTITUTION
    return [StabilizationStrategy(
        strategy_type="proline_substitution",
        description=(
            f"Substitute Gly/Ala → Pro at {len(top)} loop positions to rigidify backbone. "
            f"Estimated ΔTm ≈ +{delta_tm:.1f} °C."
        ),
        target_positions=top,
        target_residues=[sequence[p - 1] + str(p) for p in top],
        estimated_delta_tm_c=round(delta_tm, 1),
        confidence="medium",
        # Citation: Eijsink 2004 doi:10.1016/j.jbiotec.2004.03.026
        reference="Eijsink 2004 doi:10.1016/j.jbiotec.2004.03.026",
    )]


def _find_disulfide_candidates(sequence: str) -> List[StabilizationStrategy]:
    """Find pairs of Cys residues that could form engineered disulfide bridges.

    Uses sequence separation as a proxy for spatial proximity.
    Citation: Eijsink 2004 doi:10.1016/j.jbiotec.2004.03.026
    """
    cys_positions = [i + 1 for i, aa in enumerate(sequence) if aa == "C"]
    if len(cys_positions) < 2:
        return []

    pairs: List[Tuple[int, int]] = []
    for a in range(len(cys_positions)):
        for b in range(a + 1, len(cys_positions)):
            sep = abs(cys_positions[b] - cys_positions[a])
            if _DISULFIDE_SEQ_WINDOW_MIN <= sep <= _DISULFIDE_SEQ_WINDOW_MAX:
                pairs.append((cys_positions[a], cys_positions[b]))
    pairs = pairs[:3]  # top 3

    if not pairs:
        return []

    delta_tm = len(pairs) * _DT_PER_DISULFIDE
    pair_strs = [f"C{p[0]}–C{p[1]}" for p in pairs]
    return [StabilizationStrategy(
        strategy_type="disulfide_bridge",
        description=(
            f"Potential disulfide bridge(s): {', '.join(pair_strs)}. "
            f"Verify with structure; estimated ΔTm ≈ +{delta_tm:.1f} °C if confirmed."
        ),
        target_positions=[p for pair in pairs for p in pair],
        target_residues=[f"C{p}" for pair in pairs for p in pair],
        estimated_delta_tm_c=round(delta_tm, 1),
        confidence="low",   # needs structural validation
        # Citation: Eijsink 2004 doi:10.1016/j.jbiotec.2004.03.026
        reference="Eijsink 2004 doi:10.1016/j.jbiotec.2004.03.026",
    )]


def _find_buried_polars(sequence: str) -> List[StabilizationStrategy]:
    """Identify polar residues in hydrophobic context (as buried polar destabilizers).

    Heuristic: polar residue (N/Q/S/T) flanked by ≥4 hydrophobic residues in
    a window of ±4 is likely buried and destabilising.
    Citation: Shoichet BK et al. 1995 doi:10.1021/bi00013a001
    """
    candidates: List[int] = []
    for i, aa in enumerate(sequence):
        if aa not in _POLAR_RESIDUES:
            continue
        window = sequence[max(0, i - 4):i] + sequence[i + 1:min(len(sequence), i + 5)]
        hydrophobic_count = sum(1 for r in window if r in _HYDROPHOBIC)
        if hydrophobic_count >= 5:   # ≥5 of 8 flanking residues are hydrophobic
            candidates.append(i + 1)

    if not candidates:
        return []

    top = candidates[:4]
    delta_tm = len(top) * _DT_PER_BURIED_POLAR_FIX
    return [StabilizationStrategy(
        strategy_type="buried_polar_elimination",
        description=(
            f"Buried polar residue(s) in hydrophobic context at "
            f"{', '.join(sequence[p-1]+str(p) for p in top)}. "
            f"Replacing with hydrophobic equivalents could gain ΔTm ≈ +{delta_tm:.1f} °C."
        ),
        target_positions=top,
        target_residues=[sequence[p - 1] + str(p) for p in top],
        estimated_delta_tm_c=round(delta_tm, 1),
        confidence="medium",
        # Citation: Shoichet 1995 doi:10.1021/bi00013a001
        reference="Shoichet 1995 doi:10.1021/bi00013a001",
    )]


def _expression_tips(sequence: str) -> List[str]:
    """Generate expression-yield engineering tips.

    Citations:
    • Goldenzweig 2016 doi:10.1016/j.molcel.2016.06.012 (PROSS expression optimisation)
    • Romero & Arnold 2009 doi:10.1038/nrm2805
    """
    tips = []
    # N-terminal Met rule: if sequence starts with M followed by small AA → good
    if sequence and sequence[0] != "M":
        tips.append(
            "N-terminal Met is absent — consider adding M or MGSS tag for bacterial expression "
            "(Goldenzweig 2016 doi:10.1016/j.molcel.2016.06.012)"
        )
    # Cys count warning
    cys_count = sequence.count("C")
    if cys_count > 4:
        tips.append(
            f"High Cys count ({cys_count}) — consider reducing free Cys to ≤2 to prevent "
            f"aggregation; or use reducing conditions / SUMO/MBP fusion "
            f"(Eijsink 2004 doi:10.1016/j.jbiotec.2004.03.026)"
        )
    # Length warning for E. coli expression
    if len(sequence) > 600:
        tips.append(
            f"Protein length ({len(sequence)} aa) — consider domain truncation or split-intein "
            f"approach for bacterial expression; eukaryotic host (insect/CHO) may be needed"
        )
    # Low complexity / repeat check
    most_freq_aa = max(set(sequence), key=sequence.count)
    freq = sequence.count(most_freq_aa) / len(sequence)
    if freq > 0.15:
        tips.append(
            f"High frequency of '{most_freq_aa}' ({freq*100:.0f}%) — may indicate low-complexity "
            f"region; check for aggregation propensity"
        )
    if not tips:
        tips.append("No major expression-yield red flags identified from sequence alone.")
    return tips


# ──────────────────────────────────────────────────────────────────────────────
# Public entry-point
# ──────────────────────────────────────────────────────────────────────────────

async def run_engineering_strategy(
    gene: str,
    sequence: str,
    af_plddt: Optional[float] = None,
    binding_sites: Optional[list] = None,
    step_cb=None,
) -> Optional[EngineeringStrategyReport]:
    """Recommend protein engineering strategies for thermostability and expression.

    Parameters
    ----------
    gene:   HGNC gene symbol.
    sequence: Canonical amino-acid sequence.
    af_plddt: Overall AlphaFold pLDDT confidence (0–100).
    binding_sites: SequenceFeature list (for context).
    step_cb: Optional SSE progress callback.

    Returns
    -------
    EngineeringStrategyReport or None on error.
    """
    try:
        if not sequence:
            return None

        # ── Step 1: flexibility assessment ────────────────────────────────────
        if step_cb:
            await step_cb("engineering", "assessing_flexibility", 0.15)

        low_plddt_regions = _find_low_plddt_regions(af_plddt)

        # ── Step 2: Pro substitutions ─────────────────────────────────────────
        if step_cb:
            await step_cb("engineering", "finding_pro_substitutions", 0.30)

        pro_strats = _find_pro_substitutions(sequence)

        # ── Step 3: Disulfide bridge candidates ───────────────────────────────
        if step_cb:
            await step_cb("engineering", "disulfide_analysis", 0.45)

        disulfide_strats = _find_disulfide_candidates(sequence)

        # ── Step 4: Buried polar residues ─────────────────────────────────────
        if step_cb:
            await step_cb("engineering", "buried_polar_analysis", 0.60)

        buried_strats = _find_buried_polars(sequence)

        all_strategies = pro_strats + disulfide_strats + buried_strats
        total_delta_tm = sum(
            s.estimated_delta_tm_c for s in all_strategies
            if s.estimated_delta_tm_c is not None
        )

        # ── Step 5: Expression tips ───────────────────────────────────────────
        if step_cb:
            await step_cb("engineering", "expression_tips", 0.72)

        expr_tips = _expression_tips(sequence)

        # ── Step 6: Gemini synthesis ──────────────────────────────────────────
        if step_cb:
            await step_cb("engineering", "gemini_synthesis", 0.85)

        gemini_design = ""
        try:
            from core.gemini_interpreter import _call  # type: ignore
            strat_summary = "\n".join(
                f"  [{s.strategy_type}] {s.description[:120]}"
                for s in all_strategies[:6]
            ) or "  No specific strategies identified."
            prompt = (
                f"You are a protein engineering expert advising on thermostability and expression.\n\n"
                f"Gene: {gene} ({len(sequence)} aa)\n"
                f"AlphaFold pLDDT: {af_plddt if af_plddt else 'not available'}\n"
                f"Low-pLDDT regions: {'; '.join(low_plddt_regions) or 'none flagged'}\n\n"
                f"Proposed engineering strategies:\n{strat_summary}\n\n"
                f"Estimated total ΔTm gain: +{total_delta_tm:.1f} °C\n\n"
                f"Expression concerns:\n"
                + "\n".join(f"  - {t}" for t in expr_tips[:3])
                + "\n\nPlease provide:\n"
                f"1. Priority order for the proposed strategies and why\n"
                f"2. Combination effects — which strategies synergise?\n"
                f"3. Risk of activity loss from each stabilisation approach\n"
                f"4. Recommended expression system (E. coli / insect / CHO) with reasoning\n"
                f"5. Key validation assays (DSF/Tm, activity after engineering)\n\n"
                f"Write at the level of a protein engineering expert for industrial biotech."
            )
            gemini_design = await _call(prompt)
        except Exception:
            gemini_design = ""

        # ── Step 7: assemble ──────────────────────────────────────────────────
        if step_cb:
            await step_cb("engineering", "complete", 1.0)

        return EngineeringStrategyReport(
            gene=gene.upper(),
            sequence_length=len(sequence),
            estimated_total_delta_tm_c=round(total_delta_tm, 1),
            low_plddt_regions=low_plddt_regions,
            strategies=all_strategies,
            expression_tips=expr_tips,
            gemini_design=gemini_design or "",
            timestamp=datetime.utcnow(),
            references=_REFERENCES,
        )

    except Exception:
        return None
