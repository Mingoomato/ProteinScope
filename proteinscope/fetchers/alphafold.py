"""AlphaFold DB fetcher.

Fetches AlphaFold2 structure metadata and downloads PDB file locally.
No API key required — fully public.

AF3 predictions require the AlphaFold Server:
https://alphafoldserver.com
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import httpx

AF_API_BASE = "https://alphafold.ebi.ac.uk/api"
AF_FILES_BASE = "https://alphafold.ebi.ac.uk/files"
STRUCTURES_DIR = Path(__file__).parent.parent / "structures"


async def fetch_alphafold_metadata(accession: str) -> Optional[dict]:
    """Fetch AlphaFold DB prediction metadata for a UniProt accession."""
    url = f"{AF_API_BASE}/prediction/{accession}"
    async with httpx.AsyncClient(timeout=30) as client:
        r = await client.get(url)
        if r.status_code == 404:
            return None
        r.raise_for_status()
        data = r.json()
        return data[0] if isinstance(data, list) and data else None


async def download_pdb(accession: str, client: httpx.AsyncClient) -> Path:
    """Download the AlphaFold PDB file to the structures directory."""
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    pdb_url = f"{AF_FILES_BASE}/AF-{accession}-F1-model_v4.pdb"
    dest = STRUCTURES_DIR / f"{accession}.pdb"
    if dest.exists():
        return dest
    r = await client.get(pdb_url, timeout=60)
    r.raise_for_status()
    dest.write_bytes(r.content)
    return dest


async def fetch_alphafold(accession: str) -> tuple[Optional[str], Optional[float]]:
    """Fetch AlphaFold metadata and download PDB.

    Returns (pdb_url, plddt_score).
    """
    meta = await fetch_alphafold_metadata(accession)
    if not meta:
        return None, None

    pdb_url = meta.get("pdbUrl") or meta.get("cifUrl")
    plddt = meta.get("confidenceScore")

    async with httpx.AsyncClient(timeout=60) as client:
        try:
            await download_pdb(accession, client)
        except Exception:
            pass  # PDB download failing is non-fatal

    return pdb_url, plddt
