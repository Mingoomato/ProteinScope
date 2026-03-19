"""Structure comparator — TM-score and RMSD between two AlphaFold/PDB structures.

Uses biotite 1.6.0 for PDB parsing and structural superposition.
TM-score is computed using the TM-align algorithm approximation via
biotite's structural superimposition.
"""

from __future__ import annotations

import asyncio
import os
from pathlib import Path
from typing import Optional

import httpx
from pydantic import BaseModel

# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class StructuralComparison(BaseModel):
    tm_score: float           # 0–1; >0.5 = same fold
    rmsd_angstrom: float
    n_residues_compared: int
    interpretation: str       # "near-identical" / "same fold" / "different fold"


# ---------------------------------------------------------------------------
# AlphaFold download
# ---------------------------------------------------------------------------

_AF_URL = "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb"


async def download_alphafold_pdb(uniprot_id: str, save_dir: str) -> Optional[str]:
    """Download an AlphaFold v4 PDB structure for a given UniProt accession.

    Saves the file as ``{save_dir}/{uniprot_id}.pdb``.

    Args:
        uniprot_id: UniProt accession (e.g. ``P04637``).
        save_dir:   Directory to save the downloaded file.

    Returns:
        Absolute path to the saved PDB file, or None on failure.
    """
    if not uniprot_id:
        return None

    try:
        os.makedirs(save_dir, exist_ok=True)
        out_path = os.path.join(save_dir, f"{uniprot_id}.pdb")
        url = _AF_URL.format(uid=uniprot_id)

        async with httpx.AsyncClient(timeout=60) as client:
            response = await client.get(url)
            if response.status_code != 200:
                return None
            with open(out_path, "wb") as fh:
                fh.write(response.content)

        return out_path

    except Exception:
        return None


# ---------------------------------------------------------------------------
# Structural comparison
# ---------------------------------------------------------------------------

def compare_structures(
    pdb_path_a: str,
    pdb_path_b: str,
) -> Optional[StructuralComparison]:
    """Compare two PDB structures using CA-atom superimposition.

    Computes RMSD and an approximate TM-score using the TM-align formula.
    Structures are trimmed to the same number of residues before comparison.

    Args:
        pdb_path_a: Absolute path to the first PDB file.
        pdb_path_b: Absolute path to the second PDB file.

    Returns:
        StructuralComparison, or None if either file is missing or parsing fails.
    """
    if not pdb_path_a or not pdb_path_b:
        return None
    if not os.path.isfile(pdb_path_a) or not os.path.isfile(pdb_path_b):
        return None

    try:
        import biotite.structure.io.pdb as pdb_io
        import biotite.structure as struc
        import numpy as np

        # ── Parse structures ────────────────────────────────────────────────
        def _load_ca(path: str):
            pdb_file = pdb_io.PDBFile.read(path)
            structure = pdb_file.get_structure(model=1)
            # Keep only CA atoms from standard amino acids
            ca = structure[structure.atom_name == "CA"]
            return ca

        ca_a = _load_ca(pdb_path_a)
        ca_b = _load_ca(pdb_path_b)

        n_a = len(ca_a)
        n_b = len(ca_b)

        if n_a == 0 or n_b == 0:
            return None

        # Trim to the shorter length for a fair comparison
        n_compare = min(n_a, n_b)
        ca_a = ca_a[:n_compare]
        ca_b = ca_b[:n_compare]

        # ── Superimpose ──────────────────────────────────────────────────────
        fitted, _transform = struc.superimpose(ca_a, ca_b)

        # ── RMSD ─────────────────────────────────────────────────────────────
        rmsd_value = float(struc.rmsd(ca_a, fitted))

        # ── TM-score (TM-align approximation) ───────────────────────────────
        L = n_compare
        if L <= 15:
            d0 = 0.5  # avoid complex numbers for very short chains
        else:
            d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
            d0 = max(d0, 0.5)

        # Per-residue distances after superimposition
        diffs = ca_a.coord - fitted.coord          # (L, 3)
        d_sq = np.sum(diffs ** 2, axis=1)          # (L,)
        tm_sum = float(np.sum(1.0 / (1.0 + d_sq / (d0 ** 2))))
        tm_score = tm_sum / L

        # Clamp to [0, 1]
        tm_score = max(0.0, min(1.0, tm_score))

        # ── Interpretation ───────────────────────────────────────────────────
        if tm_score > 0.9:
            interpretation = "near-identical"
        elif tm_score >= 0.5:
            interpretation = "same fold"
        else:
            interpretation = "different fold"

        return StructuralComparison(
            tm_score=round(tm_score, 4),
            rmsd_angstrom=round(rmsd_value, 3),
            n_residues_compared=n_compare,
            interpretation=interpretation,
        )

    except ImportError:
        # biotite not installed
        return None
    except Exception:
        return None
