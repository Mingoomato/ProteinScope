"""KEGG fetcher — metabolic pathways, disease maps, reaction details.

KEGG REST API requires no key. Human organism prefix: hsa.
"""

from __future__ import annotations

import re
from pathlib import Path

import httpx

_BASE = "https://rest.kegg.jp"


# ---------------------------------------------------------------------------
# Pathway discovery
# ---------------------------------------------------------------------------

async def get_pathways_for_gene(ncbi_gene_id: str) -> list[str]:
    """Return KEGG pathway IDs (e.g. 'path:hsa00010') for an NCBI gene ID."""
    async with httpx.AsyncClient() as client:
        r = await client.get(f"{_BASE}/link/pathway/hsa:{ncbi_gene_id}", timeout=30)
    if r.status_code != 200 or not r.text.strip():
        return []
    ids = []
    for line in r.text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            ids.append(parts[1])
    return ids


async def get_disease_pathways_for_gene(ncbi_gene_id: str) -> list[str]:
    """Return KEGG DISEASE IDs linked to an NCBI gene ID."""
    async with httpx.AsyncClient() as client:
        r = await client.get(f"{_BASE}/link/disease/hsa:{ncbi_gene_id}", timeout=30)
    if r.status_code != 200 or not r.text.strip():
        return []
    ids = []
    for line in r.text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            ids.append(parts[1])
    return ids


# ---------------------------------------------------------------------------
# Detail fetchers
# ---------------------------------------------------------------------------

async def get_pathway_detail(pathway_id: str) -> str:
    """Return the KEGG flat-file text for a pathway entry."""
    clean_id = pathway_id.replace("path:", "")
    async with httpx.AsyncClient() as client:
        r = await client.get(f"{_BASE}/get/{clean_id}", timeout=30)
    return r.text if r.status_code == 200 else ""


async def get_reactions_for_pathway(pathway_id: str) -> list[str]:
    """Return KEGG reaction IDs linked to a pathway."""
    clean_id = pathway_id.replace("path:", "")
    async with httpx.AsyncClient() as client:
        r = await client.get(f"{_BASE}/link/reaction/{clean_id}", timeout=30)
    if r.status_code != 200 or not r.text.strip():
        return []
    ids = []
    for line in r.text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            ids.append(parts[1])
    return ids


def parse_kegg_flat(text: str) -> dict:
    """Parse a KEGG flat-file response into a dict of section → value strings."""
    result: dict[str, str] = {}
    current_key = ""
    for line in text.split("\n"):
        if line.startswith("///"):
            break
        if line and not line[0].isspace():
            parts = line.split(None, 1)
            current_key = parts[0]
            result[current_key] = parts[1] if len(parts) > 1 else ""
        elif current_key:
            result[current_key] = result[current_key] + " " + line.strip()
    return result


async def get_reaction_detail(reaction_id: str) -> dict:
    """Return parsed KEGG reaction entry (ENTRY, NAME, EQUATION, ENZYME)."""
    clean_id = reaction_id.replace("rn:", "")
    async with httpx.AsyncClient() as client:
        r = await client.get(f"{_BASE}/get/{clean_id}", timeout=30)
    if r.status_code != 200:
        return {}
    return parse_kegg_flat(r.text)


# ---------------------------------------------------------------------------
# Image downloaders
# ---------------------------------------------------------------------------

async def download_pathway_image(
    pathway_id: str, ncbi_gene_id: str, save_dir: str
) -> tuple[str | None, str | None]:
    """Download plain and gene-colored KEGG pathway PNG images.

    Returns (plain_path, colored_path). Either may be None on failure.
    """
    clean_id = pathway_id.replace("path:", "")
    plain_url = f"{_BASE}/get/{clean_id}/image"
    colored_url = (
        f"https://www.kegg.jp/kegg-bin/show_pathway?"
        f"{clean_id}&multi_query=hsa:{ncbi_gene_id}%09%23ff0000"
    )

    Path(save_dir).mkdir(parents=True, exist_ok=True)

    plain_path: str | None = None
    colored_path: str | None = None

    async with httpx.AsyncClient() as client:
        try:
            plain = await client.get(plain_url, timeout=60)
            if plain.status_code == 200 and plain.content:
                plain_path = str(Path(save_dir) / f"kegg_{clean_id}_plain.png")
                Path(plain_path).write_bytes(plain.content)
        except Exception:
            pass

        try:
            colored = await client.get(colored_url, timeout=60)
            if colored.status_code == 200 and colored.content:
                colored_path = str(Path(save_dir) / f"kegg_{clean_id}_colored.png")
                Path(colored_path).write_bytes(colored.content)
        except Exception:
            pass

    # If colored image failed, fall back to plain
    if colored_path is None:
        colored_path = plain_path

    return plain_path, colored_path


# ---------------------------------------------------------------------------
# Top-level wrapper
# ---------------------------------------------------------------------------

async def fetch_all_kegg(ncbi_gene_id: str) -> dict:
    """Return all KEGG pathway and disease data for an NCBI gene ID."""
    if not ncbi_gene_id:
        return {"pathway_ids": [], "disease_pathway_ids": [], "pathway_details": {}}

    import asyncio

    pathway_ids, disease_ids = await asyncio.gather(
        get_pathways_for_gene(ncbi_gene_id),
        get_disease_pathways_for_gene(ncbi_gene_id),
        return_exceptions=True,
    )

    if not isinstance(pathway_ids, list):
        pathway_ids = []
    if not isinstance(disease_ids, list):
        disease_ids = []

    # Fetch detail for up to 10 pathways (avoid very long PDFs)
    pathway_details = {}
    for pid in pathway_ids[:10]:
        detail_text = await get_pathway_detail(pid)
        clean_key = pid.replace("path:", "")
        pathway_details[clean_key] = parse_kegg_flat(detail_text)

    return {
        "pathway_ids": pathway_ids,
        "disease_pathway_ids": disease_ids,
        "pathway_details": pathway_details,
    }
