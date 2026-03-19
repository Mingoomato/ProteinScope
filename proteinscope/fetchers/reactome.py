"""Reactome fetcher — signaling pathways, disease cascades, pathway diagrams.

Reactome is fully open. No API key required.
"""

from __future__ import annotations

from pathlib import Path

import httpx

_BASE = "https://reactome.org/ContentService"


# ---------------------------------------------------------------------------
# Pathway discovery
# ---------------------------------------------------------------------------

async def fetch_pathways(uniprot_id: str, disease_only: bool = False) -> list[dict]:
    """Return low-level Reactome pathways for a UniProt entity.

    Pass disease_only=True to restrict to disease pathways.
    """
    params: dict = {"species": "9606"}
    if disease_only:
        params["disease"] = "true"
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/data/pathways/low/entity/{uniprot_id}",
            params=params,
            timeout=30,
        )
    if r.status_code != 200:
        return []
    try:
        return r.json()
    except Exception:
        return []


async def fetch_pathway_hierarchy(pathway_id: str) -> list[dict]:
    """Return the ancestor hierarchy (breadcrumb) for a Reactome pathway."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/data/pathway/{pathway_id}/ancestors",
            timeout=30,
        )
    if r.status_code != 200:
        return []
    try:
        return r.json()
    except Exception:
        return []


async def fetch_contained_events(pathway_id: str) -> list[dict]:
    """Return reactions/events contained within a Reactome pathway."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/data/pathway/{pathway_id}/containedEvents",
            timeout=30,
        )
    if r.status_code != 200:
        return []
    try:
        return r.json()
    except Exception:
        return []


# ---------------------------------------------------------------------------
# Image downloaders
# ---------------------------------------------------------------------------

async def download_pathway_diagram(
    pathway_id: str, uniprot_id: str, save_dir: str
) -> tuple[str | None, str | None]:
    """Download plain and protein-highlighted Reactome pathway PNGs.

    Returns (plain_path, colored_path). Either may be None on failure.
    """
    plain_url = f"{_BASE}/exporter/diagram/{pathway_id}.png?quality=8"
    colored_url = f"{_BASE}/exporter/diagram/{pathway_id}.png?quality=8&sel={uniprot_id}"

    Path(save_dir).mkdir(parents=True, exist_ok=True)

    plain_path: str | None = None
    colored_path: str | None = None

    async with httpx.AsyncClient() as client:
        try:
            plain = await client.get(plain_url, timeout=60)
            if plain.status_code == 200 and plain.content:
                plain_path = str(Path(save_dir) / f"reactome_{pathway_id}_plain.png")
                Path(plain_path).write_bytes(plain.content)
        except Exception:
            pass

        try:
            colored = await client.get(colored_url, timeout=60)
            if colored.status_code == 200 and colored.content:
                colored_path = str(Path(save_dir) / f"reactome_{pathway_id}_colored.png")
                Path(colored_path).write_bytes(colored.content)
        except Exception:
            pass

    if colored_path is None:
        colored_path = plain_path

    return plain_path, colored_path


# ---------------------------------------------------------------------------
# Top-level wrapper
# ---------------------------------------------------------------------------

async def fetch_all_reactome(uniprot_id: str) -> dict:
    """Return all Reactome pathway data (signaling + disease) for a protein."""
    import asyncio

    signaling_pathways, disease_pathways = await asyncio.gather(
        fetch_pathways(uniprot_id, disease_only=False),
        fetch_pathways(uniprot_id, disease_only=True),
        return_exceptions=True,
    )

    if not isinstance(signaling_pathways, list):
        signaling_pathways = []
    if not isinstance(disease_pathways, list):
        disease_pathways = []

    return {
        "signaling_pathways": signaling_pathways,
        "disease_pathways": disease_pathways,
    }
