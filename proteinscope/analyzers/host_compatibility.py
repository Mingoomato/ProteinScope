"""Host-Protein Expression Compatibility Matrix (P5).

Scores a protein against 7 common expression hosts based on:
- Codon Adaptation Index (CAI) from the CDS nucleotide sequence
- Disulfide bond count (from UniProt "Disulfide bond" features)
- N-glycosylation site count (from UniProt "Glycosylation" features)
- Signal peptide / transmembrane topology
- Sequence length

Zero new network calls — all data is already in the ProteinRecord.

Hosts scored:
  E.coli BL21, CHO-K1, HEK293, S.cerevisiae, Pichia pastoris, Sf9, Cell-free (WGE)
"""

from __future__ import annotations

import logging
import re
from typing import Optional, List

from pydantic import BaseModel, Field

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------


class ExpressionHostScore(BaseModel):
    host_name: str                          # e.g. "E.coli BL21(DE3)"
    host_type: str                          # "prokaryote"|"mammalian"|"yeast"|"insect"|"cell_free"
    overall_score: float                    # 0.0 – 1.0 (higher = better)
    recommendation: str                     # "Recommended"|"Acceptable"|"Caution"|"Not recommended"
    cai_score: Optional[float] = None       # Codon Adaptation Index
    disulfide_flag: str = "OK"              # "OK"|"CAUTION"|"HIGH"
    glycosylation_flag: str = "OK"          # "OK"|"CAUTION"|"NOT_SUPPORTED"
    signal_peptide_flag: str = "OK"         # "OK"|"REMOVE_REQUIRED"
    transmembrane_flag: str = "OK"          # "OK"|"CAUTION"
    size_flag: str = "OK"                   # "OK"|"LARGE"|"TOO_LARGE"
    notes: List[str] = Field(default_factory=list)


class ExpressionCompatibilityReport(BaseModel):
    gene: str
    sequence_length: int
    disulfide_bond_count: int
    glycosylation_site_count: int
    has_signal_peptide: bool
    has_transmembrane: bool
    hosts: List[ExpressionHostScore] = Field(default_factory=list)
    recommended_host: Optional[str] = None
    gemini_recommendation: Optional[str] = None  # filled by caller if needed


# ---------------------------------------------------------------------------
# CAI computation (BioPython)
# ---------------------------------------------------------------------------

# E. coli K-12 codon usage table (relative adaptiveness values, simplified)
# Source: Ikemura (1985), used as a proxy for all prokaryotes
_ECOLI_CAI_TABLE: dict[str, float] = {
    "TTT": 0.296, "TTC": 1.000, "TTA": 0.020, "TTG": 0.020,
    "CTT": 0.042, "CTC": 0.037, "CTA": 0.007, "CTG": 1.000,
    "ATT": 0.185, "ATC": 1.000, "ATA": 0.003, "ATG": 1.000,
    "GTT": 0.279, "GTC": 0.195, "GTA": 0.140, "GTG": 1.000,
    "TCT": 0.850, "TCC": 0.744, "TCA": 0.077, "TCG": 0.017,
    "CCT": 0.070, "CCC": 0.012, "CCA": 0.194, "CCG": 1.000,
    "ACT": 0.965, "ACC": 1.000, "ACA": 0.076, "ACG": 0.099,
    "GCT": 0.586, "GCC": 0.122, "GCA": 0.586, "GCG": 1.000,
    "TAT": 0.239, "TAC": 1.000, "TAA": 1.000, "TAG": 0.001,
    "CAT": 0.291, "CAC": 1.000, "CAA": 0.124, "CAG": 1.000,
    "AAT": 0.051, "AAC": 1.000, "AAA": 0.253, "AAG": 1.000,
    "GAT": 0.621, "GAC": 1.000, "GAA": 0.684, "GAG": 1.000,
    "TGT": 0.500, "TGC": 1.000, "TGA": 0.001, "TGG": 1.000,
    "CGT": 1.000, "CGC": 0.356, "CGA": 0.004, "CGG": 0.004,
    "AGT": 0.085, "AGC": 0.410, "AGA": 0.002, "AGG": 0.002,
    "GGT": 1.000, "GGC": 0.724, "GGA": 0.010, "GGG": 0.010,
}


def _compute_cai(cds: str, table: dict[str, float]) -> float:
    """Compute Codon Adaptation Index (geometric mean of relative adaptiveness)."""
    if not cds or len(cds) < 6:
        return 0.0
    cds = cds.upper().replace(" ", "").replace("\n", "")
    # Trim to reading frame (multiple of 3)
    cds = cds[: (len(cds) // 3) * 3]
    values = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3]
        if "N" in codon or len(codon) < 3:
            continue
        w = table.get(codon, 0.0)
        if w > 0:
            values.append(w)
    if not values:
        return 0.0
    import math
    log_sum = sum(math.log(v) for v in values)
    return round(math.exp(log_sum / len(values)), 4)


# ---------------------------------------------------------------------------
# Feature counting helpers
# ---------------------------------------------------------------------------


def _count_feature_type(features: list[dict], feature_type: str) -> int:
    """Count UniProt features of a given type.

    UniProt REST API v2 uses 'type' as the feature key (not 'featureType').
    """
    return sum(1 for f in features if f.get("type") == feature_type)


def _has_feature_type(features: list[dict], feature_type: str) -> bool:
    return _count_feature_type(features, feature_type) > 0


# ---------------------------------------------------------------------------
# Per-host scoring logic
# ---------------------------------------------------------------------------

def _score_host(
    host_name: str,
    host_type: str,
    disulfide_count: int,
    glyco_count: int,
    has_signal: bool,
    has_tm: bool,
    seq_len: int,
    cai: Optional[float],
    cai_table: dict[str, float],
    cds: Optional[str],
) -> ExpressionHostScore:
    notes = []
    flags = {}

    # ── Size ──
    if seq_len > 1500:
        flags["size"] = "TOO_LARGE"
        notes.append(f"Sequence length {seq_len} aa may require chaperones or truncation")
    elif seq_len > 800:
        flags["size"] = "LARGE"
        notes.append(f"Large protein ({seq_len} aa) — consider optimising expression conditions")
    else:
        flags["size"] = "OK"

    # ── Disulfide bonds ──
    if host_type == "prokaryote":
        if disulfide_count >= 3:
            flags["disulfide"] = "HIGH"
            notes.append(f"{disulfide_count} disulfide bonds — use Origami/SHuffle strains or periplasmic expression")
        elif disulfide_count >= 1:
            flags["disulfide"] = "CAUTION"
            notes.append(f"{disulfide_count} disulfide bond(s) — oxidative folding limited in E. coli cytoplasm")
        else:
            flags["disulfide"] = "OK"
    elif host_type in ("mammalian", "insect"):
        if disulfide_count >= 10:
            flags["disulfide"] = "CAUTION"
            notes.append(f"{disulfide_count} disulfide bonds — verify correct pairing")
        else:
            flags["disulfide"] = "OK"
    elif host_type == "cell_free":
        if disulfide_count >= 1:
            flags["disulfide"] = "CAUTION"
            notes.append(f"{disulfide_count} disulfide bond(s) — add DsbC + GSH/GSSG to cell-free reaction")
        else:
            flags["disulfide"] = "OK"
    else:
        flags["disulfide"] = "OK"

    # ── Glycosylation ──
    if glyco_count > 0:
        if host_type == "prokaryote":
            flags["glycosylation"] = "NOT_SUPPORTED"
            notes.append(f"{glyco_count} N-glycosylation site(s) — NOT supported in E. coli; glycoproteins require eukaryotic host")
        elif host_type == "yeast":
            flags["glycosylation"] = "CAUTION"
            notes.append(f"{glyco_count} N-glycosylation site(s) — yeast adds high-mannose glycans (may differ from human patterns)")
        elif host_type == "insect":
            flags["glycosylation"] = "CAUTION"
            notes.append(f"{glyco_count} N-glycosylation site(s) — insect cells add simple N-glycans (paucimannose), not complex human-type")
        elif host_type == "cell_free":
            flags["glycosylation"] = "NOT_SUPPORTED"
            notes.append(f"{glyco_count} N-glycosylation site(s) — cell-free systems cannot perform N-glycosylation")
        else:
            flags["glycosylation"] = "OK"  # mammalian: full glycosylation supported
    else:
        flags["glycosylation"] = "OK"

    # ── Signal peptide ──
    if has_signal and host_type == "prokaryote":
        flags["signal"] = "REMOVE_REQUIRED"
        notes.append("Signal peptide detected — remove or replace with bacterial signal sequence for periplasmic targeting")
    elif has_signal and host_type == "cell_free":
        flags["signal"] = "CAUTION"
        notes.append("Signal peptide detected — remove for cell-free expression to prevent misfolding")
    else:
        flags["signal"] = "OK"

    # ── Transmembrane ──
    if has_tm:
        if host_type in ("prokaryote", "cell_free"):
            flags["tm"] = "CAUTION"
            notes.append("Transmembrane domain — membrane protein expression may require detergent reconstitution")
        else:
            flags["tm"] = "CAUTION"
            notes.append("Transmembrane domain — verify membrane integration efficiency")
    else:
        flags["tm"] = "OK"

    # ── Overall score (0–1, weighted deductions) ──
    score = 1.0
    if flags["size"] == "TOO_LARGE":
        score -= 0.3
    elif flags["size"] == "LARGE":
        score -= 0.1
    if flags["disulfide"] == "HIGH":
        score -= 0.35
    elif flags["disulfide"] == "CAUTION":
        score -= 0.15
    if flags["glycosylation"] == "NOT_SUPPORTED":
        score -= 0.4
    elif flags["glycosylation"] == "CAUTION":
        score -= 0.15
    if flags["signal"] == "REMOVE_REQUIRED":
        score -= 0.1
    if flags["tm"] == "CAUTION":
        score -= 0.2
    # CAI bonus/malus for E. coli
    if cai is not None:
        if cai < 0.5:
            score -= 0.15
            notes.append(f"CAI={cai:.2f} — low codon adaptation; codon-optimize CDS for {host_name}")
        elif cai > 0.75:
            notes.append(f"CAI={cai:.2f} — good codon adaptation for {host_name}")

    score = max(0.0, min(1.0, score))

    if score >= 0.75:
        recommendation = "Recommended"
    elif score >= 0.5:
        recommendation = "Acceptable"
    elif score >= 0.25:
        recommendation = "Caution"
    else:
        recommendation = "Not recommended"

    return ExpressionHostScore(
        host_name=host_name,
        host_type=host_type,
        overall_score=round(score, 3),
        recommendation=recommendation,
        cai_score=cai,
        disulfide_flag=flags.get("disulfide", "OK"),
        glycosylation_flag=flags.get("glycosylation", "OK"),
        signal_peptide_flag=flags.get("signal", "OK"),
        transmembrane_flag=flags.get("tm", "OK"),
        size_flag=flags.get("size", "OK"),
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def run_expression_compatibility(
    gene: str,
    canonical_sequence: str,
    encoding_dna_cds: Optional[str],
    entry_features: list[dict],
) -> ExpressionCompatibilityReport:
    """Score protein expression compatibility across 7 standard hosts.

    This is a synchronous function (pure computation) — no network calls.
    Call it directly in query_engine.py after fetch_all() gather completes.

    Args:
        gene: Gene symbol.
        canonical_sequence: Amino acid sequence.
        encoding_dna_cds: CDS nucleotide sequence (optional — used for CAI).
        entry_features: UniProt entry["features"] list.

    Returns:
        ExpressionCompatibilityReport with per-host scoring.
    """
    if not canonical_sequence:
        return ExpressionCompatibilityReport(
            gene=gene, sequence_length=0,
            disulfide_bond_count=0, glycosylation_site_count=0,
            has_signal_peptide=False, has_transmembrane=False,
        )

    try:
        features = entry_features or []
        seq_len = len(canonical_sequence)
        disulfide_count = _count_feature_type(features, "Disulfide bond")
        glyco_count = _count_feature_type(features, "Glycosylation")
        has_signal = _has_feature_type(features, "Signal peptide")
        has_tm = _has_feature_type(features, "Transmembrane")

        # CAI for E. coli (only if CDS available)
        ecoli_cai = _compute_cai(encoding_dna_cds, _ECOLI_CAI_TABLE) if encoding_dna_cds else None

        # Define hosts
        host_specs = [
            ("E.coli BL21(DE3)", "prokaryote", ecoli_cai),
            ("CHO-K1", "mammalian", None),
            ("HEK293", "mammalian", None),
            ("S.cerevisiae", "yeast", None),
            ("Pichia pastoris", "yeast", None),
            ("Sf9 (Baculovirus)", "insect", None),
            ("Cell-free (WGE)", "cell_free", None),
        ]

        hosts = [
            _score_host(name, htype, disulfide_count, glyco_count,
                        has_signal, has_tm, seq_len, cai, _ECOLI_CAI_TABLE, encoding_dna_cds)
            for name, htype, cai in host_specs
        ]

        # Sort by score descending
        hosts.sort(key=lambda h: h.overall_score, reverse=True)
        recommended_host = hosts[0].host_name if hosts and hosts[0].overall_score >= 0.5 else None

        return ExpressionCompatibilityReport(
            gene=gene,
            sequence_length=seq_len,
            disulfide_bond_count=disulfide_count,
            glycosylation_site_count=glyco_count,
            has_signal_peptide=has_signal,
            has_transmembrane=has_tm,
            hosts=hosts,
            recommended_host=recommended_host,
        )

    except Exception as exc:
        _log.warning("Expression compatibility analysis failed for %s: %s", gene, exc)
        return ExpressionCompatibilityReport(
            gene=gene, sequence_length=len(canonical_sequence),
            disulfide_bond_count=0, glycosylation_site_count=0,
            has_signal_peptide=False, has_transmembrane=False,
        )
