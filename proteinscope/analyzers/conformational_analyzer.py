"""Conformational Ensemble Analyzer — ANM/GNM hinge detection via ProDy.

Background:
  Anisotropic Network Model (ANM): coarse-grained elastic network model.
    Each Cα atom is a node; springs connect pairs within a cutoff (typically 15 Å).
    Calculates normal modes; lowest-frequency modes describe large collective motions.
    Mean Square Fluctuations (MSF) per residue identify flexible regions.
    Local MSF minima → hinge residues that act as pivots for domain motion.

  Gaussian Network Model (GNM): isotropic variant; faster, gives B-factor predictions.
    Used to validate ANM results and identify allosteric communication paths.

Key outputs:
  - Per-residue MSF array (relative flexibility)
  - Hinge residues (local MSF minima with flanking flexibility)
  - Top 3 dominant slow modes and their functional interpretation
  - Allosteric communication: residues correlated with hinges

This analyzer:
  1. Downloads AlphaFold PDB to temp file
  2. Parses Cα atoms via ProDy (parsePDB)
  3. Runs ANM with 20 modes (Cα-only, cutoff 15 Å)
  4. Computes per-residue MSF (calcSqFlucts on first 3 slow modes)
  5. Detects hinge residues (local MSF minima, separation > 10 residues)
  6. Calls Gemini for functional interpretation
"""

from __future__ import annotations

import asyncio
import json as _json
import os
import tempfile
from typing import List, Optional

from pydantic import BaseModel, Field

from core.evidence import EvidenceGrade, DataProvenance


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

class FlexibilityMode(BaseModel):
    mode_number: int           # 1-based (mode 1 = slowest)
    eigenvalue: float          # spring constant equivalent; lower = softer
    description: str = ""      # e.g. "Large hinge motion around G719"


class ConformationalAnalysis(BaseModel):
    gene: str
    sequence_length: int
    msf_per_residue: List[float] = Field(default_factory=list)   # per-residue MSF
    hinge_residues: List[int] = Field(default_factory=list)      # 1-based positions
    high_flexibility_regions: List[str] = Field(default_factory=list)  # "aa 320-340 (IDR)"
    top_modes: List[FlexibilityMode] = Field(default_factory=list)
    prody_available: bool = False
    gemini_interpretation: str = ""
    allosteric_notes: str = ""
    # V3 Grand Consortium: ANM promoting vibration overlap (용박사, SNU/KAIST)
    promoting_vibration_modes: List[tuple] = Field(default_factory=list)
    # [(mode_index, overlap_score)] where overlap_score >= 0.5
    donor_acceptor_distance: Optional[float] = None
    provenance: Optional[DataProvenance] = None


# ---------------------------------------------------------------------------
# Hinge detection helper
# ---------------------------------------------------------------------------

def _find_hinges(msf: list[float], min_separation: int = 10) -> list[int]:
    """Find local MSF minima (hinge residues).

    A position i (0-based) is a hinge if:
    - msf[i] < msf[i-2] and msf[i] < msf[i+2]  (local minimum in 5-point window)
    - msf[i] < mean_msf * 0.7                    (below 70% of mean)
    - Separated by at least min_separation residues from any previously found hinge

    Returns 1-based positions.
    """
    if not msf or len(msf) < 5:
        return []

    n = len(msf)
    mean_msf = sum(msf) / n
    threshold = mean_msf * 0.7

    hinges: list[int] = []
    for i in range(2, n - 2):
        if (
            msf[i] < msf[i - 2]
            and msf[i] < msf[i + 2]
            and msf[i] < threshold
        ):
            # Enforce minimum separation
            if hinges and (i - (hinges[-1] - 1)) < min_separation:
                continue
            hinges.append(i + 1)  # convert to 1-based

    return hinges


# ---------------------------------------------------------------------------
# High-flexibility region detection
# ---------------------------------------------------------------------------

def _find_high_flex_regions(msf: list[float], min_len: int = 10) -> list[str]:
    """Find stretches of ≥min_len consecutive residues where msf[i] > mean*1.5.

    Returns strings like "aa 320–340".
    """
    if not msf:
        return []

    mean_msf = sum(msf) / len(msf)
    threshold = mean_msf * 1.5

    regions: list[str] = []
    start: Optional[int] = None

    for i, v in enumerate(msf):
        if v > threshold:
            if start is None:
                start = i
        else:
            if start is not None:
                length = i - start
                if length >= min_len:
                    regions.append(f"aa {start + 1}–{i}")
                start = None

    # Handle trailing region
    if start is not None:
        length = len(msf) - start
        if length >= min_len:
            regions.append(f"aa {start + 1}–{len(msf)}")

    return regions


# ---------------------------------------------------------------------------
# Gemini synthesis
# ---------------------------------------------------------------------------

async def _interpret_conformational(
    gene: str,
    n_ca: int,
    hinge_list: list[int],
    flex_regions: list[str],
    top_eigenvalue: float,
    mean_msf: float,
    max_msf: float,
) -> tuple[str, str, list[str]]:
    """Call Gemini for functional interpretation of ANM results.

    Returns (interpretation, allosteric_notes, mode_descriptions).
    """
    try:
        from core.gemini_interpreter import _call

        prompt = (
            f"Gene: {gene}\n"
            f"Sequence length (Cα atoms analyzed): {n_ca}\n"
            f"Hinge residues: {hinge_list}\n"
            f"High-flexibility regions: {flex_regions}\n"
            f"Top slow mode eigenvalue: {top_eigenvalue:.4e}\n"
            f"Mean MSF: {mean_msf:.3f}, Max MSF: {max_msf:.3f}\n\n"
            "You are a structural bioinformatics expert in normal mode analysis.\n"
            "Based on this ANM conformational data:\n"
            "1. What are the likely functional motions described by these hinge points?\n"
            "2. Are there known allosteric or activation loop regions that match the hinge positions?\n"
            "3. What are the therapeutic implications (allosteric sites, conformational epitopes)?\n\n"
            "Return ONLY raw JSON (no markdown):\n"
            '{"interpretation": "...", "allosteric_notes": "...", '
            '"mode_descriptions": ["mode 1: ...", "mode 2: ...", "mode 3: ..."]}'
        )

        raw = await _call(prompt)
        if not raw:
            return "", "", []

        cleaned = raw.strip()
        if cleaned.startswith("```"):
            cleaned = cleaned.lstrip("```json").lstrip("```").rstrip("```").strip()

        data = _json.loads(cleaned)
        interpretation = str(data.get("interpretation", "")).strip()
        allosteric_notes = str(data.get("allosteric_notes", "")).strip()
        mode_descs = data.get("mode_descriptions", [])
        if not isinstance(mode_descs, list):
            mode_descs = []
        mode_descs = [str(x).strip() for x in mode_descs if x]

        return interpretation, allosteric_notes, mode_descs

    except Exception:
        return "", "", []


# ---------------------------------------------------------------------------
# ANM Promoting Vibration Overlap (용박사, Grand Consortium V3 2026-03-20)
# ---------------------------------------------------------------------------

def calc_promoting_vibration_overlap(
    mode_vectors: list,        # list of np.ndarray, each shape (3*N,), one per mode
    ca_positions: list,        # list of (x, y, z) tuples per residue (Cα only, 0-based)
    donor_index: int,          # 0-based residue index of donor atom
    acceptor_index: int,       # 0-based residue index of acceptor atom
    threshold: float = 0.5,    # minimum overlap to include (용박사: O_k >= 0.5)
) -> list[tuple[int, float]]:
    """Compute ANM mode overlap with donor-acceptor compression vector.

    Formula (용박사, SNU/KAIST; Grand Consortium V3):
      e_DA = (R_A - R_D) / |R_A - R_D|           unit vector from donor to acceptor
      Δu_DA^k = u_A^k - u_D^k                    relative displacement in mode k
      O_k = |Δu_DA^k · e_DA| / |Δu_DA^k|         cosine overlap with DA axis

    Threshold: O_k >= 0.5 = plausible promoting mode (용박사 consensus)
    O_k >= 0.7 = strong candidate

    # Citation: Hay & Scrutton (2012) Nat Chem 4:161-168. doi:10.1038/nchem.1223
    # Citation: Hammes-Schiffer (2006) Acc Chem Res 39:93. doi:10.1021/ar040199a

    Args:
        mode_vectors: List of 3N-dimensional displacement vectors from ProDy ANM.
        ca_positions: List of (x,y,z) Cα coordinates (0-based, same order as ANM).
        donor_index: 0-based index of donor residue in ca_positions.
        acceptor_index: 0-based index of acceptor residue in ca_positions.
        threshold: Minimum overlap score to include in output.

    Returns:
        List of (mode_index_1based, overlap_score) sorted by overlap descending.
        Only modes with O_k >= threshold are returned.
    """
    try:
        import numpy as np
    except ImportError:
        return []

    if donor_index >= len(ca_positions) or acceptor_index >= len(ca_positions):
        return []
    if donor_index == acceptor_index:
        return []

    # Donor-acceptor unit vector
    r_d = np.array(ca_positions[donor_index], dtype=float)
    r_a = np.array(ca_positions[acceptor_index], dtype=float)
    da_vec = r_a - r_d
    da_norm = np.linalg.norm(da_vec)
    if da_norm < 1e-6:
        return []
    e_da = da_vec / da_norm

    results = []
    for k, mode_vec in enumerate(mode_vectors):
        try:
            vec = np.array(mode_vec, dtype=float)
            if len(vec) < (max(donor_index, acceptor_index) + 1) * 3:
                continue
            # Extract 3D displacement for donor and acceptor
            u_d = vec[donor_index * 3: donor_index * 3 + 3]
            u_a = vec[acceptor_index * 3: acceptor_index * 3 + 3]
            delta_u = u_a - u_d
            delta_norm = np.linalg.norm(delta_u)
            if delta_norm < 1e-10:
                continue
            overlap = abs(np.dot(delta_u, e_da)) / delta_norm
            if overlap >= threshold:
                results.append((k + 1, round(float(overlap), 4)))
        except Exception:
            continue

    results.sort(key=lambda x: x[1], reverse=True)
    return results


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_conformational_analysis(
    gene: str,
    af_pdb_url: str,
    step_cb=None,
) -> ConformationalAnalysis:
    """Run ANM-based conformational analysis on an AlphaFold structure.

    Args:
        gene:        Gene symbol (e.g. "EGFR", "BRCA1").
        af_pdb_url:  URL to the AlphaFold PDB file.
        step_cb:     Optional async progress callback (receives a string message).

    Returns:
        ConformationalAnalysis with MSF, hinges, flexibility regions, and
        Gemini interpretation. On any failure, returns a stub with
        prody_available=False.
    """

    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    try:
        import httpx

        # ── Step 1: Download AlphaFold PDB ───────────────────────────────────
        await _step("Downloading AlphaFold structure for conformational analysis...")

        pdb_path: Optional[str] = None
        try:
            async with httpx.AsyncClient(timeout=30) as client:
                resp = await client.get(af_pdb_url)
                resp.raise_for_status()
                with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
                    f.write(resp.content)
                    pdb_path = f.name
        except Exception as exc:
            return ConformationalAnalysis(
                gene=gene,
                sequence_length=0,
                prody_available=False,
                gemini_interpretation=f"AlphaFold PDB download failed: {exc}",
            )

        # ── Step 2: ProDy ANM ─────────────────────────────────────────────────
        await _step("Running ANM normal mode analysis (ProDy)...")

        loop = asyncio.get_event_loop()

        def _run_prody(path: str):
            import prody
            prody.confProDy(verbosity='none')
            structure = prody.parsePDB(path, subset='ca')
            if structure is None or len(structure) == 0:
                return None
            # Limit to first 700 Cα to avoid memory issues on very long sequences
            if len(structure) > 700:
                structure = structure[:700]
            # ANM
            anm = prody.ANM(f'{gene}_anm')
            anm.buildHessian(structure, cutoff=15.0)
            anm.calcModes(n_modes=min(20, len(structure) - 1))
            # Per-residue MSF from first 3 slow modes
            msf = prody.calcSqFlucts(anm[:3]).tolist()
            # Top modes info
            modes = []
            for i in range(min(3, anm.numModes())):
                m = anm[i]
                modes.append({'mode_number': i + 1, 'eigenvalue': float(m.getEigval())})
            return {'msf': msf, 'modes': modes, 'n_ca': len(structure)}

        try:
            result = await loop.run_in_executor(None, _run_prody, pdb_path)
        finally:
            try:
                os.unlink(pdb_path)
            except Exception:
                pass

        if result is None:
            return ConformationalAnalysis(
                gene=gene,
                sequence_length=0,
                prody_available=False,
                gemini_interpretation="ProDy parsed an empty structure — no Cα atoms found.",
            )

        msf: list[float] = result['msf']
        raw_modes: list[dict] = result['modes']
        n_ca: int = result['n_ca']

        # ── Step 3: Hinge residue detection ──────────────────────────────────
        await _step("Detecting hinge residues...")
        hinge_residues = _find_hinges(msf, min_separation=10)

        # ── Step 4: High-flexibility regions ─────────────────────────────────
        await _step("Identifying high-flexibility regions...")
        flex_regions = _find_high_flex_regions(msf, min_len=10)

        # Build FlexibilityMode objects (descriptions filled in by Gemini)
        top_modes = [
            FlexibilityMode(
                mode_number=m['mode_number'],
                eigenvalue=m['eigenvalue'],
            )
            for m in raw_modes
        ]

        mean_msf = sum(msf) / len(msf) if msf else 0.0
        max_msf = max(msf) if msf else 0.0
        top_ev = raw_modes[0]['eigenvalue'] if raw_modes else 0.0

        # ── Step 5: Gemini interpretation ────────────────────────────────────
        await _step("Interpreting conformational dynamics with AI...")
        interpretation, allosteric_notes, mode_descs = await _interpret_conformational(
            gene=gene,
            n_ca=n_ca,
            hinge_list=hinge_residues,
            flex_regions=flex_regions,
            top_eigenvalue=top_ev,
            mean_msf=mean_msf,
            max_msf=max_msf,
        )

        # Attach Gemini mode descriptions to FlexibilityMode objects
        for i, mode in enumerate(top_modes):
            if i < len(mode_descs):
                mode.description = mode_descs[i]

        return ConformationalAnalysis(
            gene=gene,
            sequence_length=n_ca,
            msf_per_residue=msf,
            hinge_residues=hinge_residues,
            high_flexibility_regions=flex_regions,
            top_modes=top_modes,
            prody_available=True,
            gemini_interpretation=interpretation,
            allosteric_notes=allosteric_notes,
            provenance=DataProvenance(
                source="AlphaFold EBI + ProDy ANM",
                evidence_grade=EvidenceGrade.COMPUTATIONAL,
                method="Elastic Network Model (ANM), cutoff 15 Å, 20 modes",
                scientific_caveat=(
                    "Coarse-grained Cα-only ENM; dynamics are an approximation of "
                    "true conformational space. Validate hinge predictions experimentally."
                ),
            ),
        )

    except Exception as exc:
        return ConformationalAnalysis(
            gene=gene,
            sequence_length=0,
            prody_available=False,
            gemini_interpretation=(
                f"ProDy not available — install ProDy>=2.4.0 ({exc})"
            ),
        )
