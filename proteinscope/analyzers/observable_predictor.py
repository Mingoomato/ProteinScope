"""Experimental observables prediction for NMR, FRET, and HDX-MS.

Predicts biophysical observables from sequence and AlphaFold confidence data:
  - Molecular weight
  - NMR feasibility (solution vs solid-state vs cryo-EM)
  - FRET-label candidate cysteine pairs
  - HDX-MS peptide map and fast-exchange regions

NMR size limits: Kay LE 2011 J Mol Biol doi:10.1016/j.jmb.2011.01.045
HDX-MS for protein dynamics: Englander SW 1983 Annu Rev Biochem
doi:10.1146/annurev.bi.52.070183.002521

Usage::

    from analyzers.observable_predictor import run_observable_prediction

    report = await run_observable_prediction(
        gene="BRCA1",
        sequence="MDLSALRVEEVQ...",
        af_plddt=78.4,
    )
    print(report.nmr_feasibility, report.hdx_expected_peptides)
"""

from __future__ import annotations

import re
from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel, Field

from core.evidence import DataProvenance, EvidenceGrade


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


class FRETCandidate(BaseModel):
    residue_a: int
    residue_b: int
    distance_angstrom: Optional[float] = None
    accessibility_score: float
    recommendation: str


class ObservablesPrediction(BaseModel):
    gene: str
    molecular_weight_kda: float
    fret_candidates: List[FRETCandidate] = Field(default_factory=list)
    nmr_feasibility: str          # "Solution NMR (TROSY)" | "TROSY-HSQC recommended" | "Cryo-EM or solid-state NMR preferred"
    nmr_recommendation: str
    hdx_expected_peptides: int
    hdx_fast_exchange_regions: List[str] = Field(default_factory=list)
    gemini_summary: str = ""
    timestamp: datetime


# ---------------------------------------------------------------------------
# Standard amino acid molecular weights (monoisotopic-free, average)
# Citation: Gasteiger E 2005 Protein Identification and Analysis (ExPASy)
# ---------------------------------------------------------------------------
_AA_MW: dict[str, float] = {
    "A": 89.09,   "R": 174.20,  "N": 132.12,  "D": 133.10,
    "C": 121.16,  "E": 147.13,  "Q": 146.15,  "G": 75.03,
    "H": 155.16,  "I": 131.17,  "L": 131.17,  "K": 146.19,
    "M": 149.21,  "F": 165.19,  "P": 115.13,  "S": 105.09,
    "T": 119.12,  "W": 204.23,  "Y": 181.19,  "V": 117.15,
}
_DEFAULT_AA_MW = 110.0  # fallback average Da per residue
_WATER_MW = 18.02       # subtract one water per peptide bond for full chain


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _compute_mw_kda(sequence: str) -> float:
    """Compute molecular weight in kDa from amino acid sequence.

    Uses per-residue average molecular weights from:
    Citation: Gasteiger E 2005 Protein Identification and Analysis ExPASy handbook
    """
    if not sequence:
        return 0.0
    total = sum(_AA_MW.get(aa.upper(), _DEFAULT_AA_MW) for aa in sequence)
    # Subtract water for each peptide bond (N-1 bonds → subtract (N-1) * 18.02)
    # Net: total - (len-1)*18.02
    total -= (len(sequence) - 1) * _WATER_MW
    return round(total / 1000.0, 2)


def _nmr_feasibility(mw_kda: float) -> tuple[str, str]:
    """Return (feasibility_label, recommendation) based on molecular weight.

    Citation: NMR size limits: Kay LE 2011 J Mol Biol doi:10.1016/j.jmb.2011.01.045
    """
    # Citation: NMR size limits: Kay LE 2011 J Mol Biol doi:10.1016/j.jmb.2011.01.045
    if mw_kda < 30.0:
        return (
            "Solution NMR (TROSY)",
            (
                f"At {mw_kda:.1f} kDa this protein is well within the solution NMR window. "
                "Standard TROSY or HSQC experiments on ¹⁵N/¹³C-labelled protein are feasible. "
                "Full backbone assignment should be achievable."
            ),
        )
    elif mw_kda <= 100.0:
        return (
            "TROSY-HSQC recommended",
            (
                f"At {mw_kda:.1f} kDa TROSY-HSQC with deuteration is recommended to manage "
                "transverse relaxation. Consider methyl-TROSY (ILV labelling) for side-chain "
                "dynamics. Segmental isotope labelling may be required for full assignment."
            ),
        )
    else:
        return (
            "Cryo-EM or solid-state NMR preferred",
            (
                f"At {mw_kda:.1f} kDa the protein exceeds practical solution NMR limits. "
                "Cryo-EM (sub-3 Å resolution routinely achievable) or solid-state NMR (MAS) "
                "are preferred structural methods. NMR is still useful for local dynamics "
                "studies on isolated domains or sparse labelling schemes."
            ),
        )


def _find_cys_positions(sequence: str) -> list[int]:
    """Return 1-based positions of all cysteine residues."""
    return [i + 1 for i, aa in enumerate(sequence.upper()) if aa == "C"]


def _estimate_inter_residue_distance(pos_a: int, pos_b: int) -> float:
    """Estimate Cα–Cα distance from sequence separation.

    For a disordered/random-coil region: ~3.8 Å / 3 per residue separation.
    This is a lower bound; structured regions may be shorter (helix ~1.5 Å/res
    along helix axis) or longer (extended β-strand ~3.5 Å/res).

    Citation: Average Cα-Cα distance in disordered regions:
    Flory PJ 1969 Statistical Mechanics of Chain Molecules
    Using 3.8/3 ≈ 1.27 Å per residue separation as disordered estimate.
    """
    sep = abs(pos_a - pos_b)
    # Citation: Flory PJ 1969 Statistical Mechanics of Chain Molecules
    distance = sep * (3.8 / 3.0)
    return round(distance, 1)


def _fret_candidates(
    sequence: str,
    af_plddt: Optional[float],
    top_n: int = 5,
) -> list[FRETCandidate]:
    """Find top FRET-label candidate Cys pairs in 30–80 Å range.

    Uses sequence-separation estimate for Cα–Cα distance.
    Accessibility_score = 1.0 if af_plddt > 70 (structured, surface-accessible)
    else 0.5 (disordered or low-confidence region).

    Citation: FRET distance range: Stryer L 1978 Annu Rev Biochem
    doi:10.1146/annurev.bi.47.070178.000525
    """
    cys_positions = _find_cys_positions(sequence)
    if len(cys_positions) < 2:
        return []

    # Citation: FRET distance range 30-80 Å: Stryer L 1978 Annu Rev Biochem doi:10.1146/annurev.bi.47.070178.000525
    fret_min_a = 30.0
    fret_max_a = 80.0

    # Accessibility: use global pLDDT as proxy; per-residue not available here
    accessibility = 1.0 if (af_plddt is not None and af_plddt > 70.0) else 0.5

    candidates: list[tuple[int, int, float]] = []
    for i in range(len(cys_positions)):
        for j in range(i + 1, len(cys_positions)):
            pos_a = cys_positions[i]
            pos_b = cys_positions[j]
            dist = _estimate_inter_residue_distance(pos_a, pos_b)
            if fret_min_a <= dist <= fret_max_a:
                candidates.append((pos_a, pos_b, dist))

    # Sort by distance closest to 50 Å (near R₀ for typical FRET pairs)
    candidates.sort(key=lambda x: abs(x[2] - 50.0))

    result: list[FRETCandidate] = []
    for pos_a, pos_b, dist in candidates[:top_n]:
        result.append(
            FRETCandidate(
                residue_a=pos_a,
                residue_b=pos_b,
                distance_angstrom=dist,
                accessibility_score=accessibility,
                recommendation=(
                    f"Cys{pos_a}–Cys{pos_b}: estimated distance {dist:.1f} Å "
                    f"(FRET window 30–80 Å). "
                    f"Accessibility score={accessibility:.1f}. "
                    "Recommend maleimide-reactive FRET dye pair (e.g. Cy3/Cy5)."
                ),
            )
        )
    return result


def _mock_trypsin_digest(sequence: str) -> list[str]:
    """Mock trypsin digest: cleave after K or R, but not if followed by P.

    Returns list of peptide strings.

    Citation: Trypsin specificity: Olsen JV 2004 Mol Cell Proteomics
    doi:10.1074/mcp.T400003-MCP200
    """
    if not sequence:
        return []

    # Citation: Trypsin specificity (K/R not before P):
    # Olsen JV 2004 Mol Cell Proteomics doi:10.1074/mcp.T400003-MCP200
    peptides: list[str] = []
    current: list[str] = []
    seq_upper = sequence.upper()

    for i, aa in enumerate(seq_upper):
        current.append(aa)
        if aa in ("K", "R"):
            # Do not cleave if next residue is P (missed cleavage rule)
            if i + 1 < len(seq_upper) and seq_upper[i + 1] == "P":
                continue
            peptides.append("".join(current))
            current = []

    if current:
        peptides.append("".join(current))

    return [p for p in peptides if p]


def _fast_exchange_regions(sequence: str, window: int = 5) -> list[str]:
    """Identify fast-exchange regions: stretches of >= `window` consecutive residues
    with > 50% polar/charged residues (S/T/D/E/N/Q).

    Fast-exchanging backbone amides in polar/disordered regions are characteristic
    of high solvent accessibility.

    Citation: HDX-MS for protein dynamics: Englander SW 1983 Annu Rev Biochem
    doi:10.1146/annurev.bi.52.070183.002521
    """
    # Citation: HDX-MS for protein dynamics: Englander SW 1983 Annu Rev Biochem doi:10.1146/annurev.bi.52.070183.002521
    fast_residues = set("STDNEQ")
    seq_upper = sequence.upper()
    n = len(seq_upper)
    regions: list[str] = []

    i = 0
    while i < n:
        # Find start of a potential fast-exchange stretch
        stretch_start = i
        stretch: list[int] = []
        j = i
        while j < n:
            # Check local window
            if seq_upper[j] in fast_residues:
                stretch.append(j)
                j += 1
            elif j + 1 < n and seq_upper[j + 1] in fast_residues:
                # Allow one non-fast-exchange residue within window
                stretch.append(j)
                j += 1
            else:
                break

        if len(stretch) >= window:
            start_1based = stretch[0] + 1
            end_1based = stretch[-1] + 1
            # Check > 50% polar in this stretch
            subseq = seq_upper[stretch[0]: stretch[-1] + 1]
            polar_frac = sum(1 for aa in subseq if aa in fast_residues) / len(subseq)
            if polar_frac > 0.5:
                regions.append(f"residues {start_1based}–{end_1based} ({subseq})")
            i = stretch[-1] + 1
        else:
            i += 1

    return regions


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


async def run_observable_prediction(
    gene: str,
    sequence: str,
    af_plddt: Optional[float] = None,
    step_cb=None,
) -> Optional[ObservablesPrediction]:
    """Predict biophysical observables for NMR, FRET, and HDX-MS.

    Args:
        gene:      Gene symbol (e.g. "BRCA1").
        sequence:  Full amino acid sequence string.
        af_plddt:  Mean AlphaFold pLDDT score (0–100); used for accessibility scoring.
        step_cb:   Optional async progress callback (receives str).

    Returns:
        ObservablesPrediction on success, None on unrecoverable failure.
    """

    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    try:
        if not sequence:
            return None

        # ── Step 1: Compute molecular weight ────────────────────────────────
        await _step("[1/4] Computing molecular weight...")

        mw_kda = _compute_mw_kda(sequence)

        # ── Step 2: Predict NMR feasibility ─────────────────────────────────
        await _step("[2/4] Predicting NMR feasibility...")

        # Citation: NMR size limits: Kay LE 2011 J Mol Biol doi:10.1016/j.jmb.2011.01.045
        nmr_feasibility, nmr_recommendation = _nmr_feasibility(mw_kda)

        # ── Step 3: Find FRET-label candidates ──────────────────────────────
        await _step("[3/4] Finding FRET-label candidates...")

        fret_candidates = _fret_candidates(sequence, af_plddt)

        # ── Step 4: HDX-MS peptide prediction ───────────────────────────────
        await _step("[4/4] HDX-MS peptide prediction...")

        # Citation: HDX-MS for protein dynamics: Englander SW 1983 Annu Rev Biochem doi:10.1146/annurev.bi.52.070183.002521
        peptides = _mock_trypsin_digest(sequence)
        hdx_expected_peptides = len(peptides)
        fast_regions = _fast_exchange_regions(sequence)

        # Gemini summary
        gemini_text = ""
        try:
            from core.gemini_interpreter import _call

            fret_summary = (
                f"{len(fret_candidates)} candidate Cys pair(s) in FRET range (30–80 Å)"
                if fret_candidates
                else "No Cys-pair FRET candidates found (no Cys residues or none in range)"
            )
            fast_summary = (
                f"{len(fast_regions)} fast-exchange region(s): {'; '.join(fast_regions[:3])}"
                if fast_regions
                else "No fast-exchange regions detected"
            )

            prompt = (
                f"Gene: {gene}\n"
                f"Molecular weight: {mw_kda:.1f} kDa\n"
                f"NMR feasibility: {nmr_feasibility}\n"
                f"FRET: {fret_summary}\n"
                f"HDX-MS: {hdx_expected_peptides} tryptic peptides expected; {fast_summary}\n\n"
                "You are a structural biophysics expert. Summarize:\n"
                "1. Which biophysical technique (NMR, FRET, HDX-MS, cryo-EM) is most"
                " informative for this protein and why?\n"
                "2. What key conformational questions can these techniques answer?\n"
                "3. Any critical experimental caveats?\n\n"
                "Return a concise 3-4 sentence summary."
            )

            gemini_text = await _call(prompt) or ""
        except Exception:
            gemini_text = ""

        return ObservablesPrediction(
            gene=gene,
            molecular_weight_kda=mw_kda,
            fret_candidates=fret_candidates,
            nmr_feasibility=nmr_feasibility,
            nmr_recommendation=nmr_recommendation,
            hdx_expected_peptides=hdx_expected_peptides,
            hdx_fast_exchange_regions=fast_regions,
            gemini_summary=gemini_text.strip(),
            timestamp=datetime.utcnow(),
        )

    except Exception:
        return None
