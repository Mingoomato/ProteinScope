"""PharmGKB fetcher — drug interactions, PGx variants, FDA labels, pathway maps.

Gold-standard pharmacogenomics database. Requires a free API key:
register at https://www.pharmgkb.org/page/contacts
"""

from __future__ import annotations

import os
from pathlib import Path

import httpx
from dotenv import load_dotenv

load_dotenv()

BASE = "https://api.pharmgkb.org/v1/data"
_KEY = os.getenv("PHARMGKB_API_KEY", "")
_HEADERS = {"Authorization": f"Token {_KEY}"} if _KEY else {}


async def get_pharmgkb_id(gene_symbol: str) -> str | None:
    """Resolve an HGNC gene symbol to a PharmGKB gene ID (e.g. 'PA35')."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{BASE}/gene",
            params={"symbol": gene_symbol, "view": "min"},
            headers=_HEADERS,
            timeout=30,
        )
        if r.status_code != 200:
            return None
        results = r.json().get("data", [])
        return results[0]["id"] if results else None


async def fetch_drug_interactions(pharmgkb_gene_id: str) -> list[dict]:
    """Return all chemicals (drugs) related to the given PharmGKB gene ID."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{BASE}/gene/{pharmgkb_gene_id}/relatedChemicals",
            headers=_HEADERS,
            timeout=30,
        )
        if r.status_code != 200:
            return []
        return r.json().get("data", [])


async def fetch_pgx_variants(pharmgkb_gene_id: str) -> list[dict]:
    """Return clinical pharmacogenomic annotations for a gene."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{BASE}/clinicalAnnotation",
            params={"gene": pharmgkb_gene_id, "view": "max"},
            headers=_HEADERS,
            timeout=30,
        )
        if r.status_code != 200:
            return []
        return r.json().get("data", [])


async def fetch_pharmgkb_pathways(pharmgkb_gene_id: str) -> list[dict]:
    """Return PharmGKB pathway entries that involve this gene."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{BASE}/pathway",
            params={"relatedGene": pharmgkb_gene_id},
            headers=_HEADERS,
            timeout=30,
        )
        if r.status_code != 200:
            return []
        return r.json().get("data", [])


async def download_pathway_image(pathway_id: str, save_dir: str) -> str | None:
    """Download a PharmGKB pathway PNG. Returns the local file path or None on failure."""
    url = f"https://api.pharmgkb.org/v1/data/pathway/{pathway_id}/download"
    async with httpx.AsyncClient() as client:
        r = await client.get(
            url, params={"format": "png"}, headers=_HEADERS, timeout=60
        )
    if r.status_code != 200 or not r.content:
        return None
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    path = str(Path(save_dir) / f"pharmgkb_{pathway_id}.png")
    with open(path, "wb") as f:
        f.write(r.content)
    return path


async def fetch_fda_labels(pharmgkb_gene_id: str) -> list[dict]:
    """Return FDA drug label annotations that mention this gene."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{BASE}/label",
            params={"relatedGene": pharmgkb_gene_id},
            headers=_HEADERS,
            timeout=30,
        )
        if r.status_code != 200:
            return []
        return r.json().get("data", [])


async def fetch_all_pharmgkb(gene_symbol: str) -> dict:
    """Top-level convenience wrapper. Returns a dict with all PharmGKB data."""
    pgkb_id = await get_pharmgkb_id(gene_symbol)
    if not pgkb_id:
        return {
            "pharmgkb_id": None,
            "drug_interactions": [],
            "pgx_variants": [],
            "pathways": [],
            "fda_labels": [],
        }

    import asyncio

    drug_interactions, pgx_variants, pathways, fda_labels = await asyncio.gather(
        fetch_drug_interactions(pgkb_id),
        fetch_pgx_variants(pgkb_id),
        fetch_pharmgkb_pathways(pgkb_id),
        fetch_fda_labels(pgkb_id),
        return_exceptions=True,
    )

    return {
        "pharmgkb_id": pgkb_id,
        "drug_interactions": drug_interactions if isinstance(drug_interactions, list) else [],
        "pgx_variants": pgx_variants if isinstance(pgx_variants, list) else [],
        "pathways": pathways if isinstance(pathways, list) else [],
        "fda_labels": fda_labels if isinstance(fda_labels, list) else [],
    }
