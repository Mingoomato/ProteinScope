"""OMIM API fetcher for disease and deficiency data.

Requires a free academic API key stored in OMIM_API_KEY env var.
Register at: https://www.omim.org/api
"""

from __future__ import annotations

import os
from typing import Optional

import httpx

OMIM_BASE = "https://api.omim.org/api"


def _api_key() -> Optional[str]:
    return os.getenv("OMIM_API_KEY")


async def search_omim(gene_name: str, client: httpx.AsyncClient) -> list[dict]:
    """Search OMIM entries for a gene name."""
    key = _api_key()
    if not key:
        return []
    params = {
        "search": gene_name,
        "include": "all",
        "format": "json",
        "start": 0,
        "limit": 10,
        "apiKey": key,
    }
    r = await client.get(f"{OMIM_BASE}/entry/search", params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    return data.get("omim", {}).get("searchResponse", {}).get("entryList", [])


async def fetch_omim_entry(mim_number: str, client: httpx.AsyncClient) -> dict:
    """Fetch a specific OMIM entry by MIM number."""
    key = _api_key()
    if not key:
        return {}
    params = {
        "mimNumber": mim_number,
        "include": "text",
        "format": "json",
        "apiKey": key,
    }
    r = await client.get(f"{OMIM_BASE}/entry", params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    entries = data.get("omim", {}).get("entryList", [])
    return entries[0].get("entry", {}) if entries else {}


def parse_diseases(entry: dict) -> list[str]:
    """Extract disease names from an OMIM entry dict."""
    diseases = []
    title = entry.get("titles", {}).get("preferredTitle", "")
    if title:
        diseases.append(title)
    text_sections = entry.get("textSectionList", [])
    for section in text_sections:
        ts = section.get("textSection", {})
        if ts.get("textSectionName") in ("clinicalFeatures", "description"):
            content = ts.get("textSectionContent", "")
            if content and len(content) > 20:
                # Just record the section title as a pointer
                diseases.append(f"[See OMIM {entry.get('mimNumber','')} for {ts.get('textSectionName','')}]")
    return diseases


async def fetch_omim_diseases(gene_name: str) -> tuple[list[str], list[str]]:
    """Return (disease_names, deficiency_diseases) for a gene."""
    if not _api_key():
        return [], []
    async with httpx.AsyncClient() as client:
        results = await search_omim(gene_name, client)
        diseases, deficiencies = [], []
        for entry_wrapper in results[:5]:
            entry_meta = entry_wrapper.get("entry", {})
            mim = str(entry_meta.get("mimNumber", ""))
            if not mim:
                continue
            entry = await fetch_omim_entry(mim, client)
            names = parse_diseases(entry)
            diseases.extend(names)
            # heuristic: if title contains "deficiency" flag as deficiency disease
            for n in names:
                if "deficiency" in n.lower() or "defect" in n.lower():
                    deficiencies.append(n)
        return diseases, deficiencies
