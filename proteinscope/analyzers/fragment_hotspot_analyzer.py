"""Fragment-Based Drug Discovery (FBDD) Hotspot Analyzer.

Identifies fragment-sized binding pockets (100-300 Å³) on AlphaFold structures
using fpocket, scores them for FBDD suitability via a Fragsite-proxy score,
checks active site overlap, and synthesises an FBDD strategy via Gemini.

Falls back gracefully to a Gemini-only analysis if fpocket is not installed
or the PDB download fails.

async def run_fragment_hotspot_analysis(gene, af_pdb_url, binding_sites, active_sites, step_cb)
    -> Optional[FragmentHotspotMap]
"""

from __future__ import annotations

import asyncio
import json as _json
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from pydantic import BaseModel, Field

from core.evidence import EvidenceGrade, DataProvenance


# ---------------------------------------------------------------------------
# Fragment pocket size constants for FBDD
# Citation: Fragment pocket size for FBDD: Kozakov D 2015 J Chem Inf Model
#           doi:10.1021/acs.jcim.5b00523
# ---------------------------------------------------------------------------

_FRAGMENT_VOLUME_MIN = 100.0    # Å³ — lower bound for fragment-sized pocket
_FRAGMENT_VOLUME_MAX = 300.0    # Å³ — upper bound for fragment-sized pocket
_OPTIMAL_FRAGMENT_VOLUME = 200.0  # Å³ — ideal target for Fragsite score formula


def _fpocket_binary() -> str:
    """Return fpocket binary path from env var or default 'fpocket'."""
    return os.environ.get("FPOCKET_BINARY", "fpocket")


def _is_binary_available(binary: str) -> bool:
    """Check if the fpocket binary is on PATH or is a valid absolute path."""
    return shutil.which(binary) is not None


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

class FragmentHotspot(BaseModel):
    pocket_id: str
    fragsite_score: float
    volume_angstrom3: float
    center_x: float
    center_y: float
    center_z: float
    residues: List[str]
    overlaps_active_site: bool
    druggability_class: str     # "FBDD-suitable" | "Too small" | "Too large"
    provenance: Optional[DataProvenance] = None


class FragmentHotspotMap(BaseModel):
    gene: str
    hotspots: List[FragmentHotspot]
    top_fragment_pocket: Optional[FragmentHotspot] = None
    fpocket_available: bool = True
    gemini_fbdd_strategy: str = ""
    timestamp: datetime
    provenance: Optional[DataProvenance] = None


# ---------------------------------------------------------------------------
# PDB download
# ---------------------------------------------------------------------------

async def _download_pdb(url: str, dest_path: str) -> bool:
    """Download AlphaFold PDB to dest_path using httpx.

    Returns True on success. Does not raise.
    """
    try:
        import httpx
        async with httpx.AsyncClient(timeout=60) as client:
            r = await client.get(url, follow_redirects=True)
            r.raise_for_status()
            Path(dest_path).write_bytes(r.content)
        return True
    except Exception:
        return False


# ---------------------------------------------------------------------------
# fpocket output parsing
# ---------------------------------------------------------------------------

def _parse_info_txt(info_path: Path) -> list[dict]:
    """Parse fpocket *_info.txt into a list of pocket property dicts.

    Each dict: {pocket_id, score, druggability_score, volume_A3, hydrophobicity}
    """
    text = info_path.read_text(encoding="utf-8", errors="ignore")
    blocks = re.split(r"Pocket\s+(\d+)\s*:", text)

    pockets: list[dict] = []
    i = 1
    while i + 1 < len(blocks):
        pocket_id = int(blocks[i].strip())
        content = blocks[i + 1]
        i += 2

        def _extract_float(pattern: str, text: str) -> Optional[float]:
            m = re.search(pattern, text)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    return None
            return None

        score = _extract_float(r"Score\s*:\s*([\d.+-]+)", content)
        druggability_score = _extract_float(r"Druggability Score\s*:\s*([\d.+-]+)", content)
        volume = _extract_float(r"Volume\s*:\s*([\d.+-]+)", content)
        hydrophobicity = _extract_float(r"Mean local hydrophobicity density\s*:\s*([\d.+-]+)", content)
        # Center of mass coordinates (if present)
        cx = _extract_float(r"x_barycenter\s*:\s*([\d.+-]+)", content)
        cy = _extract_float(r"y_barycenter\s*:\s*([\d.+-]+)", content)
        cz = _extract_float(r"z_barycenter\s*:\s*([\d.+-]+)", content)

        pockets.append({
            "pocket_id": pocket_id,
            "score": score if score is not None else 0.0,
            "druggability_score": druggability_score if druggability_score is not None else 0.0,
            "volume_A3": volume if volume is not None else 0.0,
            "hydrophobicity": hydrophobicity if hydrophobicity is not None else 0.0,
            "center_x": cx if cx is not None else 0.0,
            "center_y": cy if cy is not None else 0.0,
            "center_z": cz if cz is not None else 0.0,
            "residues": [],
        })

    return pockets


def _parse_pocket_residues(pocket_pdb_path: Path) -> list[str]:
    """Extract unique residue identifiers from a pocket *_atm.pdb file.

    Returns list like ["ALA42", "GLY55", ...] (residue name + sequence number).
    """
    residues: list[str] = []
    seen: set[str] = set()
    try:
        text = pocket_pdb_path.read_text(encoding="utf-8", errors="ignore")
        for line in text.splitlines():
            if not line.startswith("ATOM"):
                continue
            try:
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                key = f"{res_name}{res_seq}"
                if key and key not in seen:
                    seen.add(key)
                    residues.append(key)
            except IndexError:
                continue
    except Exception:
        pass
    return residues


# ---------------------------------------------------------------------------
# Scoring and classification
# ---------------------------------------------------------------------------

def _fragsite_score(druggability_score: float, volume: float) -> float:
    """Compute Fragsite score proxy for a pocket.

    Formula: score = druggability_score * (1 - abs(volume - 200) / 200)
    Higher is better; peaks at 200 Å³ (ideal fragment pocket volume).
    Clamped to [0, 1].

    Citation: Fragment pocket size for FBDD: Kozakov D 2015 J Chem Inf Model
    doi:10.1021/acs.jcim.5b00523
    """
    if volume <= 0:
        return 0.0
    shape_factor = 1.0 - abs(volume - _OPTIMAL_FRAGMENT_VOLUME) / _OPTIMAL_FRAGMENT_VOLUME
    shape_factor = max(0.0, shape_factor)
    raw = float(druggability_score) * shape_factor
    return max(0.0, min(1.0, raw))


def _druggability_class(volume: float) -> str:
    """Classify pocket size for FBDD suitability.

    Citation: Fragment pocket size for FBDD: Kozakov D 2015 J Chem Inf Model
    doi:10.1021/acs.jcim.5b00523
    """
    if volume < _FRAGMENT_VOLUME_MIN:
        return "Too small"
    if volume > _FRAGMENT_VOLUME_MAX:
        return "Too large"
    return "FBDD-suitable"


# ---------------------------------------------------------------------------
# Active site overlap
# ---------------------------------------------------------------------------

_NUM_RE = re.compile(r"(\d+)$")


def _extract_residue_numbers(sites: list) -> set[int]:
    """Extract integer residue positions from binding_sites / active_sites lists.

    Handles UniProt feature dicts, SequenceFeature objects, and plain int/str.
    """
    positions: set[int] = set()
    if not sites:
        return positions

    for item in sites:
        try:
            if isinstance(item, int):
                positions.add(item)
                continue
            if isinstance(item, str):
                m = _NUM_RE.search(item)
                if m:
                    positions.add(int(m.group(1)))
                continue
            # UniProt feature dict
            if isinstance(item, dict):
                loc = item.get("location", {})
                start = loc.get("start", {})
                end = loc.get("end", {})
                try:
                    s_val = int(start.get("value", 0)) if isinstance(start, dict) else int(start)
                    e_val = int(end.get("value", s_val)) if isinstance(end, dict) else int(end)
                    for pos in range(s_val, e_val + 1):
                        if pos > 0:
                            positions.add(pos)
                except (TypeError, ValueError):
                    pass
                continue
            # SequenceFeature-like object
            for attr in ("start", "position", "begin"):
                val = getattr(item, attr, None)
                if val is not None:
                    try:
                        positions.add(int(val))
                    except (TypeError, ValueError):
                        pass
                    break
        except Exception:
            continue

    return positions


def _check_active_site_overlap(pocket_residues: list[str], active_positions: set[int]) -> bool:
    """Return True if any pocket residue position is in the active site position set."""
    if not active_positions or not pocket_residues:
        return False
    for res in pocket_residues:
        m = _NUM_RE.search(res)
        if m:
            try:
                if int(m.group(1)) in active_positions:
                    return True
            except ValueError:
                continue
    return False


# ---------------------------------------------------------------------------
# Gemini synthesis
# ---------------------------------------------------------------------------

async def _synthesize_fbdd_strategy(
    gene: str,
    hotspots: List[FragmentHotspot],
    fpocket_available: bool,
    binding_sites: list,
    active_sites: list,
) -> str:
    """Call Gemini for FBDD strategy synthesis.

    When fpocket is not available, Gemini reasons from gene name and known
    binding site annotations alone. Returns strategy string; empty on failure.
    """
    try:
        from core.gemini_interpreter import _call

        if fpocket_available and hotspots:
            top_spots = sorted(hotspots, key=lambda h: h.fragsite_score, reverse=True)[:3]
            hotspot_lines = "\n".join(
                f"  - {h.pocket_id}: Fragsite={h.fragsite_score:.3f}, "
                f"vol={h.volume_angstrom3:.0f} Å³, class={h.druggability_class}, "
                f"active_site_overlap={h.overlaps_active_site}"
                for h in top_spots
            )
            pocket_context = f"Top fragment hotspots:\n{hotspot_lines}"
        else:
            bs_note = f"{len(binding_sites)} binding site annotation(s)" if binding_sites else "no binding site data"
            as_note = f"{len(active_sites)} active site annotation(s)" if active_sites else "no active site data"
            pocket_context = (
                f"fpocket not available — reasoning from gene name and annotations only. "
                f"Known annotations: {bs_note}, {as_note}."
            )

        prompt = (
            f"Gene: {gene}\n"
            f"{pocket_context}\n\n"
            "You are an expert in fragment-based drug discovery (FBDD). "
            "Based on the above fragment hotspot data:\n"
            "1. Which fragment hotspot(s) are most promising for FBDD and why?\n"
            "2. What fragment screening strategy is recommended "
            "(biophysical methods: SPR, STD-NMR, thermal shift, X-ray fragment soaking)?\n"
            "3. What are the key risks for fragment elaboration at this target "
            "(pocket flexibility, water displacement, linker geometry, MW growth)?\n"
            "4. Are there opportunities for fragment-to-lead via "
            "structure-guided growing, merging, or linking?\n\n"
            "Return ONLY raw JSON (no markdown):\n"
            '{"best_hotspot": "...", "screening_strategy": "...", '
            '"elaboration_approach": "...", "key_risks": ["...", "..."], '
            '"synthesis": "2-3 sentences"}'
        )

        raw = await _call(prompt)
        if not raw:
            return ""

        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)

        parts = []
        synthesis = str(data.get("synthesis", "")).strip()
        if synthesis:
            parts.append(synthesis)
        screening = str(data.get("screening_strategy", "")).strip()
        if screening:
            parts.append(f"Screening: {screening}")
        elaboration = str(data.get("elaboration_approach", "")).strip()
        if elaboration:
            parts.append(f"Elaboration: {elaboration}")
        risks = data.get("key_risks", [])
        if isinstance(risks, list) and risks:
            parts.append("Risks: " + "; ".join(str(r) for r in risks if r))

        return " | ".join(parts) if parts else ""

    except Exception:
        return ""


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_fragment_hotspot_analysis(
    gene: str,
    af_pdb_url: str,
    binding_sites: list,
    active_sites: list,
    step_cb=None,
) -> Optional[FragmentHotspotMap]:
    """Run fragment hotspot analysis for a target protein.

    Args:
        gene:          Gene symbol (e.g. "EGFR").
        af_pdb_url:    URL to AlphaFold PDB file (EBI AlphaFold API).
        binding_sites: List of UniProt binding site features or SequenceFeature objects.
        active_sites:  List of UniProt active site features or SequenceFeature objects.
        step_cb:       Optional async progress callback (receives str message).

    Returns:
        FragmentHotspotMap — None only on catastrophic unrecoverable failure.
    """
    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    tmpdir: Optional[str] = None
    fpocket_out_dir: Optional[Path] = None

    try:
        # ── Step 1: Download AlphaFold PDB ────────────────────────────────────
        await _step("[1/5] Downloading AlphaFold PDB...")

        tmpdir = tempfile.mkdtemp(prefix="frag_hotspot_")
        pdb_path = os.path.join(tmpdir, "structure.pdb")

        ok = await _download_pdb(af_pdb_url, pdb_path)
        if not ok or not os.path.exists(pdb_path):
            gemini_strategy = await _synthesize_fbdd_strategy(
                gene, [], False, binding_sites, active_sites
            )
            return FragmentHotspotMap(
                gene=gene,
                hotspots=[],
                fpocket_available=False,
                gemini_fbdd_strategy="PDB download failed" if not gemini_strategy else gemini_strategy,
                timestamp=datetime.utcnow(),
                provenance=DataProvenance(
                    source="ProteinScope Fragment Hotspot Analyzer (Gemini-only)",
                    evidence_grade=EvidenceGrade.AI_GENERATED,
                    scientific_caveat="PDB download failed; no structural hotspot analysis performed.",
                    method="Gemini LLM synthesis only",
                ),
            )

        # ── Step 2: Run fpocket ───────────────────────────────────────────────
        await _step("[2/5] Running fpocket for fragment-sized pockets...")

        binary = _fpocket_binary()
        fpocket_avail = _is_binary_available(binary)

        pockets_raw: list[dict] = []

        if not fpocket_avail:
            # fpocket not installed — fall through to Gemini-only analysis
            await _step("[3/5] Scoring fragment hotspots...")
            await _step("[4/5] Checking active site overlap...")
            await _step("[5/5] Gemini FBDD synthesis (fpocket unavailable)...")

            gemini_strategy = await _synthesize_fbdd_strategy(
                gene, [], False, binding_sites, active_sites
            )
            return FragmentHotspotMap(
                gene=gene,
                hotspots=[],
                fpocket_available=False,
                gemini_fbdd_strategy=gemini_strategy,
                timestamp=datetime.utcnow(),
                provenance=DataProvenance(
                    source="ProteinScope Fragment Hotspot Analyzer (Gemini-only)",
                    evidence_grade=EvidenceGrade.AI_GENERATED,
                    scientific_caveat=(
                        "fpocket binary not available; fragment hotspot geometry not computed. "
                        "Install fpocket for structural analysis."
                    ),
                    method="Gemini LLM synthesis only",
                ),
            )

        # fpocket is available — run it
        try:
            proc = await asyncio.create_subprocess_exec(
                binary, "-f", pdb_path,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=tmpdir,
            )
            try:
                stdout, stderr = await asyncio.wait_for(
                    proc.communicate(), timeout=120
                )
            except asyncio.TimeoutError:
                # fpocket timed out — degrade gracefully
                proc.kill()
                fpocket_avail = False

        except FileNotFoundError:
            fpocket_avail = False
        except Exception:
            fpocket_avail = False

        if fpocket_avail and proc.returncode != 0:
            fpocket_avail = False

        if fpocket_avail:
            # Locate fpocket output directory
            pdb_stem = Path(pdb_path).stem   # "structure"
            out_dir = Path(tmpdir) / f"{pdb_stem}_out"
            if not out_dir.exists():
                candidates = list(Path(tmpdir).glob("*_out"))
                out_dir = candidates[0] if candidates else None

            fpocket_out_dir = out_dir

            if out_dir and out_dir.exists():
                info_files = list(out_dir.glob("*_info.txt"))
                if info_files:
                    try:
                        pockets_raw = _parse_info_txt(info_files[0])
                    except Exception:
                        pockets_raw = []

                    # Enrich pocket residues from pocket PDB files
                    pockets_dir = out_dir / "pockets"
                    if pockets_dir.exists() and pockets_raw:
                        for pocket in pockets_raw:
                            pid = pocket["pocket_id"]
                            atm_pdb = pockets_dir / f"pocket{pid}_atm.pdb"
                            if atm_pdb.exists():
                                pocket["residues"] = _parse_pocket_residues(atm_pdb)

        # Filter for fragment-sized pockets (100-300 Å³)
        # Citation: Fragment pocket size for FBDD: Kozakov D 2015 J Chem Inf Model
        # doi:10.1021/acs.jcim.5b00523
        fragment_pockets = [
            p for p in pockets_raw
            if _FRAGMENT_VOLUME_MIN <= float(p.get("volume_A3", 0.0)) <= _FRAGMENT_VOLUME_MAX
        ]

        # ── Step 3: Score fragment hotspots ──────────────────────────────────
        await _step("[3/5] Scoring fragment hotspots...")

        # Citation: Fragment pocket size for FBDD: Kozakov D 2015 J Chem Inf Model
        # doi:10.1021/acs.jcim.5b00523
        # score = druggability_score * (1 - abs(volume - 200) / 200)
        scored_hotspots: list[dict] = []
        for p in fragment_pockets:
            vol = float(p.get("volume_A3", 0.0))
            drug_score = float(p.get("druggability_score", 0.0))
            fs = _fragsite_score(drug_score, vol)
            cls = _druggability_class(vol)
            scored_hotspots.append({
                "pocket_id": f"pocket{p['pocket_id']}",
                "fragsite_score": fs,
                "volume_angstrom3": vol,
                "center_x": float(p.get("center_x", 0.0)),
                "center_y": float(p.get("center_y", 0.0)),
                "center_z": float(p.get("center_z", 0.0)),
                "residues": p.get("residues", []),
                "druggability_class": cls,
            })

        # ── Step 4: Check active site overlap ────────────────────────────────
        await _step("[4/5] Checking active site overlap...")

        # Extract canonical active site positions from binding_sites and active_sites
        active_positions = _extract_residue_numbers(binding_sites or [])
        active_positions |= _extract_residue_numbers(active_sites or [])

        hotspot_models: List[FragmentHotspot] = []
        for h in scored_hotspots:
            overlap = _check_active_site_overlap(h["residues"], active_positions)
            hotspot_models.append(
                FragmentHotspot(
                    pocket_id=h["pocket_id"],
                    fragsite_score=h["fragsite_score"],
                    volume_angstrom3=h["volume_angstrom3"],
                    center_x=h["center_x"],
                    center_y=h["center_y"],
                    center_z=h["center_z"],
                    residues=h["residues"],
                    overlaps_active_site=overlap,
                    druggability_class=h["druggability_class"],
                    provenance=DataProvenance(
                        source="fpocket + AlphaFold EBI v4",
                        evidence_grade=EvidenceGrade.COMPUTATIONAL,
                        scientific_caveat=(
                            "Fragment hotspot geometry derived from AlphaFold backbone; "
                            "Fragsite score is a computational proxy — experimental "
                            "fragment screening required for validation."
                        ),
                        method=(
                            "fpocket geometric pocket detection; "
                            "Fragsite score = druggability_score × (1 − |vol−200|/200)"
                        ),
                    ),
                )
            )

        # Sort by Fragsite score descending
        hotspot_models.sort(key=lambda h: h.fragsite_score, reverse=True)
        top_hotspot = hotspot_models[0] if hotspot_models else None

        # ── Step 5: Gemini synthesis ──────────────────────────────────────────
        await _step("[5/5] Gemini FBDD synthesis...")

        gemini_strategy = await _synthesize_fbdd_strategy(
            gene, hotspot_models, fpocket_avail, binding_sites, active_sites
        )

        return FragmentHotspotMap(
            gene=gene,
            hotspots=hotspot_models,
            top_fragment_pocket=top_hotspot,
            fpocket_available=fpocket_avail,
            gemini_fbdd_strategy=gemini_strategy,
            timestamp=datetime.utcnow(),
            provenance=DataProvenance(
                source="fpocket + AlphaFold EBI v4 + Gemini synthesis",
                evidence_grade=EvidenceGrade.COMPUTATIONAL,
                scientific_caveat=(
                    "Fragment hotspot geometry from AlphaFold; fpocket scores are "
                    "predictive. Experimental fragment screening (X-ray soaking, "
                    "SPR, STD-NMR) required for confirmation."
                ),
                method="fpocket pocket detection filtered to fragment-sized pockets (100-300 Å³)",
            ),
        )

    except Exception:
        return None

    finally:
        # Always clean up tempfile and fpocket output directory
        if tmpdir:
            try:
                shutil.rmtree(tmpdir, ignore_errors=True)
            except Exception:
                pass
