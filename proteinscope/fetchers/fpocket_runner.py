"""fpocket binding pocket runner.

Uses the fpocket binary (system install). Reads FPOCKET_BINARY env var,
defaults to "fpocket". Downloads AlphaFold PDB to tempfile, runs fpocket,
parses _info.txt output. Returns [] gracefully if binary absent or fails.

async def run_fpocket(af_pdb_url: str) -> list[dict]
# Returns: [{pocket_id, score, druggability_score, volume_A3, hydrophobicity, residues[]}]
"""

from __future__ import annotations

import asyncio
import os
import re
import shutil
import tempfile
from pathlib import Path
from typing import Optional


def _fpocket_binary() -> str:
    return os.environ.get("FPOCKET_BINARY", "fpocket")


def _is_binary_available(binary: str) -> bool:
    """Check if the fpocket binary is on PATH or is a valid absolute path."""
    return shutil.which(binary) is not None


async def _download_pdb(url: str, dest_path: str) -> bool:
    """Download a PDB file from url to dest_path. Returns True on success."""
    try:
        import httpx
        async with httpx.AsyncClient(timeout=60) as client:
            r = await client.get(url, follow_redirects=True)
            r.raise_for_status()
            Path(dest_path).write_bytes(r.content)
        return True
    except Exception:
        return False


def _parse_info_txt(info_path: Path) -> list[dict]:
    """Parse fpocket *_info.txt into a list of pocket dicts.

    Each dict: {pocket_id, score, druggability_score, volume_A3, hydrophobicity}
    Residues are populated separately from pocket PDB files.
    """
    text = info_path.read_text(encoding="utf-8", errors="ignore")

    # Split into pocket blocks — each block starts with "Pocket N :"
    blocks = re.split(r"Pocket\s+(\d+)\s*:", text)
    # blocks[0] = header before first pocket
    # blocks[1], blocks[2] = (pocket_num, content), (pocket_num, content), ...

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

        pockets.append({
            "pocket_id": pocket_id,
            "score": score if score is not None else 0.0,
            "druggability_score": druggability_score if druggability_score is not None else 0.0,
            "volume_A3": volume if volume is not None else 0.0,
            "hydrophobicity": hydrophobicity if hydrophobicity is not None else 0.0,
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
            # PDB ATOM format:
            # Columns 18-20: residue name (3-char)
            # Columns 23-26: residue sequence number
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


async def run_fpocket(af_pdb_url: str) -> list[dict]:
    """Run fpocket on an AlphaFold PDB and return parsed pocket data.

    Args:
        af_pdb_url: URL to an AlphaFold PDB file (e.g. from EBI AlphaFold API).

    Returns:
        List of pocket dicts, each with keys:
            pocket_id (int), score (float), druggability_score (float),
            volume_A3 (float), hydrophobicity (float), residues (list[str])
        Returns [] on any failure (binary absent, download fail, parse error).
    """
    binary = _fpocket_binary()

    # Graceful degradation — return [] if fpocket binary is not installed
    if not _is_binary_available(binary):
        return []

    tmpdir = tempfile.mkdtemp(prefix="fpocket_")
    try:
        pdb_path = os.path.join(tmpdir, "structure.pdb")

        # Step 1: Download PDB
        ok = await _download_pdb(af_pdb_url, pdb_path)
        if not ok or not os.path.exists(pdb_path):
            return []

        # Step 2: Run fpocket
        try:
            proc = await asyncio.create_subprocess_exec(
                binary, "-f", pdb_path,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=tmpdir,
            )
            stdout, stderr = await asyncio.wait_for(proc.communicate(), timeout=120)
        except asyncio.TimeoutError:
            return []
        except FileNotFoundError:
            return []
        except Exception:
            return []

        if proc.returncode != 0:
            return []

        # Step 3: Locate output directory — fpocket creates structure_out/
        out_dir = Path(tmpdir) / "structure_out"
        if not out_dir.exists():
            # Sometimes the name follows the pdb filename without extension
            candidates = list(Path(tmpdir).glob("*_out"))
            if candidates:
                out_dir = candidates[0]
            else:
                return []

        # Step 4: Parse info file
        info_files = list(out_dir.glob("*_info.txt"))
        if not info_files:
            return []
        info_path = info_files[0]

        try:
            pockets = _parse_info_txt(info_path)
        except Exception:
            return []

        if not pockets:
            return []

        # Step 5: Enrich each pocket with residue list from pocket PDB files
        pockets_dir = out_dir / "pockets"
        if pockets_dir.exists():
            for pocket in pockets:
                pid = pocket["pocket_id"]
                atm_pdb = pockets_dir / f"pocket{pid}_atm.pdb"
                if atm_pdb.exists():
                    pocket["residues"] = _parse_pocket_residues(atm_pdb)

        return pockets

    except Exception:
        return []
    finally:
        # Always clean up temp directory
        try:
            shutil.rmtree(tmpdir, ignore_errors=True)
        except Exception:
            pass
