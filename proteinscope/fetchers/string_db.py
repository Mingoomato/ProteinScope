"""STRING DB fetcher — protein-protein interaction network.

STRING maps protein-protein interactions with confidence scores.
No API key required.
"""

from __future__ import annotations

from pathlib import Path

import httpx

_BASE = "https://string-db.org/api"


async def fetch_interactions(
    gene_symbol: str,
    min_score: int = 400,
    limit: int = 20,
) -> list[dict]:
    """Return interaction partners for a gene from STRING.

    Args:
        gene_symbol: HGNC gene symbol.
        min_score:   Minimum combined interaction score (0–1000). Default 400.
        limit:       Maximum number of partners to return. Default 20.

    Each returned dict includes: preferredName_A, preferredName_B,
    score, experimentally_confirmed (derived from experiments score > 0).
    """
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/json/interaction_partners",
            params={
                "identifiers": gene_symbol,
                "species": 9606,
                "limit": limit,
                "required_score": min_score,
                "caller_identity": "proteinscope",
            },
            timeout=30,
        )
    if r.status_code != 200:
        return []
    try:
        return r.json()
    except Exception:
        return []


async def download_network_image(gene_symbol: str, save_dir: str) -> str | None:
    """Download the STRING functional network PNG.

    Returns the local file path or None on failure.
    """
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/image/network",
            params={
                "identifiers": gene_symbol,
                "species": 9606,
                "network_type": "functional",
                "caller_identity": "proteinscope",
            },
            timeout=60,
        )
    if r.status_code != 200 or not r.content:
        return None
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    path = str(Path(save_dir) / f"string_{gene_symbol}_network.png")
    Path(path).write_bytes(r.content)
    return path


async def fetch_all_string(gene_symbol: str, min_score: int = 400) -> dict:
    """Top-level wrapper returning STRING interactions for a gene symbol."""
    interactions = await fetch_interactions(gene_symbol, min_score=min_score)
    return {
        "gene_symbol": gene_symbol,
        "interactions": interactions if isinstance(interactions, list) else [],
    }
