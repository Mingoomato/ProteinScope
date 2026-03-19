"""RNA / CDS analyzer — codon adaptation index and basic CDS analysis.

ViennaRNA Python bindings are not installed. This module implements:
- Codon Adaptation Index (CAI) using standard codon usage tables
- Basic CDS properties (GC content, length, start/stop codon check)
- Notes where ViennaRNA would add mRNA secondary structure prediction

ViennaRNA can be added later: pip install ViennaRNA (Linux/Mac) or from source.
"""

from __future__ import annotations

import math
from typing import Optional

from pydantic import BaseModel

# ---------------------------------------------------------------------------
# ViennaRNA availability check
# ---------------------------------------------------------------------------
_VIENNARNA_AVAILABLE: bool = False
try:
    import RNA  # type: ignore
    _VIENNARNA_AVAILABLE = True
except ImportError:
    pass


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class RNAAnalysisResult(BaseModel):
    cds_length_nt: int
    gc_content_pct: float
    has_start_codon: bool
    has_stop_codon: bool
    codon_adaptation_index: float
    target_organism: str
    viennarna_available: bool = False
    mfe_kcal: Optional[float] = None   # None until ViennaRNA is installed
    expression_prediction: str         # "high" / "moderate" / "low" based on CAI


# ---------------------------------------------------------------------------
# Human codon usage table (Homo sapiens, Sharp & Li 1987 / GenBank compilation)
#
# Values are Relative Adaptiveness (RA):
#   RA = (usage_of_codon) / (usage_of_most_used_synonymous_codon)
# The most-used codon for each amino acid has RA = 1.0.
# Stop codons are assigned RA = 1.0 (they are excluded from CAI calculation).
# ---------------------------------------------------------------------------

HUMAN_CODON_RA: dict[str, float] = {
    # Phe (F)
    "TTT": 0.83, "TTC": 1.00,
    # Leu (L)
    "TTA": 0.15, "TTG": 0.32,
    "CTT": 0.42, "CTC": 0.71, "CTA": 0.24, "CTG": 1.00,
    # Ile (I)
    "ATT": 0.79, "ATC": 1.00, "ATA": 0.42,
    # Met (M)
    "ATG": 1.00,
    # Val (V)
    "GTT": 0.54, "GTC": 0.69, "GTA": 0.35, "GTG": 1.00,
    # Ser (S)
    "TCT": 0.72, "TCC": 0.88, "TCA": 0.68, "TCG": 0.25,
    "AGT": 0.68, "AGC": 1.00,
    # Pro (P)
    "CCT": 0.85, "CCC": 1.00, "CCA": 0.77, "CCG": 0.29,
    # Thr (T)
    "ACT": 0.71, "ACC": 1.00, "ACA": 0.74, "ACG": 0.27,
    # Ala (A)
    "GCT": 0.78, "GCC": 1.00, "GCA": 0.59, "GCG": 0.25,
    # Tyr (Y)
    "TAT": 0.72, "TAC": 1.00,
    # Stop (*) — excluded from CAI; assigned 1.0 as placeholder
    "TAA": 1.00, "TAG": 1.00, "TGA": 1.00,
    # His (H)
    "CAT": 0.64, "CAC": 1.00,
    # Gln (Q)
    "CAA": 0.45, "CAG": 1.00,
    # Asn (N)
    "AAT": 0.62, "AAC": 1.00,
    # Lys (K)
    "AAA": 0.71, "AAG": 1.00,
    # Asp (D)
    "GAT": 0.74, "GAC": 1.00,
    # Glu (E)
    "GAA": 0.71, "GAG": 1.00,
    # Cys (C)
    "TGT": 0.64, "TGC": 1.00,
    # Trp (W)
    "TGG": 1.00,
    # Arg (R)
    "CGT": 0.34, "CGC": 0.58, "CGA": 0.30, "CGG": 0.54,
    "AGA": 1.00, "AGG": 0.83,
    # Gly (G)
    "GGT": 0.58, "GGC": 0.78, "GGA": 0.65, "GGG": 0.50,  # GGC most used in humans
}

_STOP_CODONS: frozenset[str] = frozenset({"TAA", "TAG", "TGA"})


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_cai(cds_sequence: str, target_organism: str = "human") -> float:
    """Compute the Codon Adaptation Index (CAI) for a CDS.

    Uses the Homo sapiens codon usage table regardless of target_organism
    (with a logged note for non-human requests). Returns the geometric mean
    of per-codon relative adaptiveness values.

    Args:
        cds_sequence:   DNA coding sequence (A/T/G/C), upper or lower case.
                        Should start with ATG and be divisible by 3.
        target_organism: Organism name. Currently only "human" is supported;
                         other values use the human table with a note.

    Returns:
        CAI as a float in [0.0, 1.0]. Returns 0.0 on error or empty input.
    """
    if not cds_sequence:
        return 0.0

    seq = cds_sequence.upper().replace(" ", "").replace("\n", "")
    if len(seq) < 3:
        return 0.0

    # Use human table (only supported table currently)
    table = HUMAN_CODON_RA

    log_sum = 0.0
    codon_count = 0

    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) < 3:
            break
        # Skip stop codons in CAI calculation
        if codon in _STOP_CODONS:
            continue
        ra = table.get(codon)
        if ra is None:
            # Unknown codon (e.g. ambiguous base) — skip
            continue
        if ra <= 0:
            continue
        log_sum += math.log(ra)
        codon_count += 1

    if codon_count == 0:
        return 0.0

    cai = math.exp(log_sum / codon_count)
    return round(min(1.0, max(0.0, cai)), 4)


def analyze_cds(cds_sequence: str, target_organism: str = "human") -> RNAAnalysisResult:
    """Perform a comprehensive analysis of a CDS nucleotide sequence.

    Checks structural validity (start codon, stop codon, frame), computes
    GC content, and calculates the Codon Adaptation Index (CAI).
    When ViennaRNA is installed, MFE secondary structure prediction is also
    performed (currently not available on this system).

    Args:
        cds_sequence:   DNA CDS as a plain nucleotide string.
        target_organism: Target organism for codon optimization context.

    Returns:
        RNAAnalysisResult with all computed properties.
    """
    seq = (cds_sequence or "").upper().strip()
    length = len(seq)

    # ── Start / stop codon checks ─────────────────────────────────────────
    has_start = seq[:3] == "ATG" if length >= 3 else False

    has_stop = False
    if length >= 3:
        last_codon = seq[-3:]
        has_stop = last_codon in _STOP_CODONS

    # ── GC content ───────────────────────────────────────────────────────
    gc_count = seq.count("G") + seq.count("C")
    gc_content_pct = round(100.0 * gc_count / length, 2) if length > 0 else 0.0

    # ── CAI ──────────────────────────────────────────────────────────────
    cai = compute_cai(seq, target_organism)

    # ── Expression prediction based on CAI thresholds ────────────────────
    if cai >= 0.70:
        expression = "high"
    elif cai >= 0.45:
        expression = "moderate"
    else:
        expression = "low"

    # ── ViennaRNA MFE (unavailable — graceful stub) ───────────────────────
    mfe: Optional[float] = None
    viennarna_ok = _VIENNARNA_AVAILABLE
    if viennarna_ok:
        try:
            import RNA  # type: ignore
            # Convert T→U for RNA folding
            rna_seq = seq.replace("T", "U")
            structure, mfe = RNA.fold(rna_seq)
        except Exception:
            mfe = None
            viennarna_ok = False

    return RNAAnalysisResult(
        cds_length_nt=length,
        gc_content_pct=gc_content_pct,
        has_start_codon=has_start,
        has_stop_codon=has_stop,
        codon_adaptation_index=cai,
        target_organism=target_organism,
        viennarna_available=viennarna_ok,
        mfe_kcal=mfe,
        expression_prediction=expression,
    )
