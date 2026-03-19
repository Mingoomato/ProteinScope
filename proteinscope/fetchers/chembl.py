"""ChEMBL fetcher — drug targets, mechanism of action, binding data.

ChEMBL provides drug binding mechanisms and IC50/Ki data. No API key required.
"""

from __future__ import annotations

import httpx

_BASE = "https://www.ebi.ac.uk/chembl/api/data"


async def get_target_chembl_id(uniprot_id: str) -> str | None:
    """Resolve a UniProt accession to a ChEMBL target ID."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/target",
            params={
                "target_components__accession": uniprot_id,
                "format": "json",
                "limit": 1,
            },
            timeout=30,
        )
        if r.status_code != 200:
            return None
        targets = r.json().get("targets", [])
        if not targets:
            return None
        return targets[0].get("target_chembl_id")


async def fetch_mechanisms(target_chembl_id: str) -> list[dict]:
    """Return mechanism-of-action entries for a ChEMBL target."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/mechanism",
            params={"target_chembl_id": target_chembl_id, "format": "json", "limit": 100},
            timeout=30,
        )
        if r.status_code != 200:
            return []
        return r.json().get("mechanisms", [])


async def fetch_activities(target_chembl_id: str, limit: int = 50) -> list[dict]:
    """Return drug activity records (IC50, Ki, Kd) for a ChEMBL target."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/activity",
            params={
                "target_chembl_id": target_chembl_id,
                "format": "json",
                "limit": limit,
                "standard_type__in": "IC50,Ki,Kd,EC50",
            },
            timeout=30,
        )
        if r.status_code != 200:
            return []
        return r.json().get("activities", [])


async def fetch_molecule(chembl_id: str) -> dict:
    """Return molecule metadata (preferred name, clinical phase, etc.)."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/molecule/{chembl_id}",
            params={"format": "json"},
            timeout=30,
        )
        if r.status_code != 200:
            return {}
        return r.json()


async def fetch_all_chembl(uniprot_id: str) -> dict:
    """Top-level wrapper. Returns mechanisms and key activity data for a protein."""
    target_id = await get_target_chembl_id(uniprot_id)
    if not target_id:
        return {"target_chembl_id": None, "mechanisms": [], "activities": []}

    import asyncio

    mechanisms, activities = await asyncio.gather(
        fetch_mechanisms(target_id),
        fetch_activities(target_id),
        return_exceptions=True,
    )

    return {
        "target_chembl_id": target_id,
        "mechanisms": mechanisms if isinstance(mechanisms, list) else [],
        "activities": activities if isinstance(activities, list) else [],
    }
