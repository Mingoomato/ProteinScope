"""Binding Pocket Druggability Suite (P6).

Integrates fpocket pocket detection with UniProt active/binding site data and
Gemini synthesis to produce a DruggabilityReport for a target protein.

async def run_druggability_analysis(gene, af_pdb_url, uniprot_entry, step_cb=None) -> DruggabilityReport
"""

from __future__ import annotations

import json as _json
from typing import List, Optional

from pydantic import BaseModel, Field

from core.evidence import EvidenceGrade, DataProvenance


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

class BindingPocket(BaseModel):
    pocket_id: int
    score: float
    druggability_score: float
    volume_A3: float
    hydrophobicity: float
    residues: List[str]
    druggability_class: str          # "Druggable" | "Tractable" | "Challenging"
    active_site_overlap: bool = False
    notes: List[str] = Field(default_factory=list)


class DruggabilityReport(BaseModel):
    gene: str
    pockets_found: int
    top_pocket: Optional[BindingPocket] = None
    all_pockets: List[BindingPocket] = Field(default_factory=list)
    fpocket_available: bool = False
    druggability_class: str = "Unknown"      # overall class for the protein
    gemini_synthesis: str = ""
    inhibition_strategies: List[str] = Field(default_factory=list)
    provenance: Optional[DataProvenance] = None


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------

def _classify_pocket(druggability_score: float) -> str:
    """Classify a pocket based on its fpocket druggability score."""
    if druggability_score >= 0.5:
        return "Druggable"
    if druggability_score >= 0.2:
        return "Tractable"
    return "Challenging"


def _overall_class(pockets: list[BindingPocket]) -> str:
    """Derive protein-level druggability class from pocket list."""
    if not pockets:
        return "Unknown"
    classes = [p.druggability_class for p in pockets]
    if "Druggable" in classes:
        return "Druggable"
    if "Tractable" in classes:
        return "Tractable"
    return "Challenging"


# ---------------------------------------------------------------------------
# Active site overlap detection
# ---------------------------------------------------------------------------

def _extract_active_site_positions(uniprot_entry: dict) -> set[int]:
    """Return residue sequence positions annotated as Active site or Binding site in UniProt."""
    positions: set[int] = set()
    for feature in uniprot_entry.get("features", []):
        ft = feature.get("type", "").lower()
        if ft in ("active site", "binding site"):
            loc = feature.get("location", {})
            # UniProt location can be a single position or a range
            start = loc.get("start", {})
            end = loc.get("end", {})
            try:
                s_val = int(start.get("value", 0))
                e_val = int(end.get("value", s_val))
                for pos in range(s_val, e_val + 1):
                    if pos > 0:
                        positions.add(pos)
            except (TypeError, ValueError):
                pass
    return positions


def _check_active_site_overlap(pocket_residues: list[str], active_positions: set[int]) -> bool:
    """Check whether any pocket residue position matches a UniProt active/binding site position.

    pocket_residues format: ["ALA42", "GLY55", ...] (residue name + seq number).
    """
    if not active_positions or not pocket_residues:
        return False
    import re
    _num_re = re.compile(r"(\d+)$")
    for res in pocket_residues:
        m = _num_re.search(res)
        if m:
            try:
                pos = int(m.group(1))
                if pos in active_positions:
                    return True
            except ValueError:
                continue
    return False


# ---------------------------------------------------------------------------
# Gemini synthesis
# ---------------------------------------------------------------------------

async def _synthesize_druggability(
    gene: str,
    pockets: list[BindingPocket],
    report: "DruggabilityReport",
) -> "DruggabilityReport":
    """Call Gemini for druggability assessment and inhibition strategies."""
    try:
        from core.gemini_interpreter import _call

        top = report.top_pocket
        if top:
            top_line = (
                f"Top pocket: druggability={top.druggability_score:.2f}, "
                f"volume={top.volume_A3:.0f} Å³, "
                f"{len(top.residues)} residues, "
                f"class={top.druggability_class}"
            )
            active_site_str = "yes" if top.active_site_overlap else "no"
        else:
            top_line = "No pockets detected"
            active_site_str = "no"

        prompt = (
            f"Gene: {gene}\n"
            f"Pockets found: {len(pockets)}\n"
            f"{top_line}\n"
            f"Active site overlap: {active_site_str}\n\n"
            "You are a structural bioinformatics expert. Based on this fpocket druggability data:\n"
            "1. What is the most promising binding site and why?\n"
            "2. What small-molecule inhibition strategies are most feasible?\n"
            "3. What are the key caveats (AlphaFold backbone reliability, pocket flexibility)?\n\n"
            "Return ONLY raw JSON (no markdown):\n"
            '{"druggability_summary": "...", "inhibition_strategies": ["...", "..."], '
            '"synthesis": "2-3 sentences"}'
        )

        raw = await _call(prompt)
        if not raw:
            return report

        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)

        report.gemini_synthesis = str(data.get("synthesis", "") or data.get("druggability_summary", "")).strip()
        strategies = data.get("inhibition_strategies", [])
        if isinstance(strategies, list):
            report.inhibition_strategies = [str(s).strip() for s in strategies if s]

    except Exception:
        pass

    return report


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_druggability_analysis(
    gene: str,
    af_pdb_url: str,
    uniprot_entry: dict,
    step_cb=None,
) -> DruggabilityReport:
    """Run binding pocket druggability analysis for a protein.

    Args:
        gene:          Gene symbol (e.g. "EGFR").
        af_pdb_url:    URL to AlphaFold PDB file.
        uniprot_entry: Raw UniProt entry dict (for active site annotation).
        step_cb:       Optional async progress callback (receives str message).

    Returns:
        DruggabilityReport — always returns, never raises. Empty report on failure.
    """
    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    report = DruggabilityReport(
        gene=gene,
        pockets_found=0,
        provenance=DataProvenance(
            source="fpocket + AlphaFold EBI v4",
            evidence_grade=EvidenceGrade.COMPUTATIONAL,
            scientific_caveat=(
                "Pocket geometry derived from AlphaFold backbone; "
                "druggability scores are predictive and should be validated experimentally."
            ),
            method="fpocket geometric pocket detection",
        ),
    )

    try:
        # ── Step 1: Download structure (handled inside fpocket_runner) ──────
        await _step("[1/5] Downloading AlphaFold structure...")

        # ── Step 2: Run fpocket ──────────────────────────────────────────────
        await _step("[2/5] Running fpocket pocket detection...")

        from fetchers.fpocket_runner import run_fpocket, _is_binary_available, _fpocket_binary
        pockets_raw = await run_fpocket(af_pdb_url)

        binary = _fpocket_binary()
        report.fpocket_available = _is_binary_available(binary)

        if not pockets_raw:
            # fpocket not available or returned no pockets — Gemini-only path
            await _step("[3/5] Classifying pockets...")
            await _step("[4/5] Checking active site overlap...")
            await _step("[5/5] Synthesizing druggability assessment...")
            report = await _synthesize_druggability(gene, [], report)
            return report

        # ── Step 3: Classify pockets ─────────────────────────────────────────
        await _step("[3/5] Classifying pockets...")

        binding_pockets: list[BindingPocket] = []
        for raw in pockets_raw:
            drug_score = float(raw.get("druggability_score", 0.0))
            cls = _classify_pocket(drug_score)
            pocket = BindingPocket(
                pocket_id=int(raw.get("pocket_id", 0)),
                score=float(raw.get("score", 0.0)),
                druggability_score=drug_score,
                volume_A3=float(raw.get("volume_A3", 0.0)),
                hydrophobicity=float(raw.get("hydrophobicity", 0.0)),
                residues=raw.get("residues", []),
                druggability_class=cls,
                notes=[],
            )
            binding_pockets.append(pocket)

        # Sort by druggability score descending
        binding_pockets.sort(key=lambda p: p.druggability_score, reverse=True)

        # ── Step 4: Active site overlap ──────────────────────────────────────
        await _step("[4/5] Checking active site overlap...")

        active_positions = _extract_active_site_positions(uniprot_entry)
        for pocket in binding_pockets:
            pocket.active_site_overlap = _check_active_site_overlap(
                pocket.residues, active_positions
            )
            if pocket.active_site_overlap:
                pocket.notes.append("Overlaps UniProt active/binding site annotation")

        report.all_pockets = binding_pockets
        report.pockets_found = len(binding_pockets)
        report.top_pocket = binding_pockets[0] if binding_pockets else None
        report.druggability_class = _overall_class(binding_pockets)

        # ── Step 5: Gemini synthesis ──────────────────────────────────────────
        await _step("[5/5] Synthesizing druggability assessment...")
        report = await _synthesize_druggability(gene, binding_pockets, report)

    except Exception:
        # Graceful degradation — return whatever report has been built so far
        pass

    return report
