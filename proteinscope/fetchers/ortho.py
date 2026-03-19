"""OrthoDB fetcher for cross-species ortholog lookup.

Uses the OrthoDB REST API to find orthologs of a given UniProt accession
across taxonomic groups.
No API key required.
"""

from __future__ import annotations

from typing import Optional

import httpx

ORTHODB_BASE = "https://data.orthodb.org/v12"


async def search_orthodb(gene_name: str, organism_taxid: str = "9606") -> Optional[str]:
    """Search OrthoDB for a gene and return the ortholog group ID."""
    async with httpx.AsyncClient(timeout=30) as client:
        params = {
            "query": gene_name,
            "ncbi_tax_id": organism_taxid,
            "limit": 1,
        }
        r = await client.get(f"{ORTHODB_BASE}/search", params=params)
        if r.status_code != 200:
            return None
        data = r.json()
        results = data.get("data", [])
        if results:
            return results[0]  # ortholog group ID
    return None


async def fetch_orthologs(group_id: str, level_taxid: str = "33208") -> list[dict]:
    """Fetch all members of an ortholog group at a given taxonomic level."""
    async with httpx.AsyncClient(timeout=30) as client:
        params = {"id": group_id, "limit": 50}
        r = await client.get(f"{ORTHODB_BASE}/orthologs", params=params)
        if r.status_code != 200:
            return []
        data = r.json()
        members = data.get("data", [])
        orthologs = []
        for member in members:
            org_name = member.get("organism", {}).get("name", "")
            gene_id = member.get("gene_id", {}).get("id", "")
            xrefs = member.get("xrefs", {})
            uniprot_id = ""
            for xref in xrefs.get("UniProt", []):
                uniprot_id = xref
                break
            if org_name:
                orthologs.append({
                    "organism": org_name,
                    "gene_id": gene_id,
                    "uniprot_id": uniprot_id,
                })
        return orthologs


async def fetch_cross_species_entries(gene_name: str, taxid: str = "9606") -> list[dict]:
    """Top-level: return ortholog list for a gene."""
    group_id = await search_orthodb(gene_name, taxid)
    if not group_id:
        return []
    return await fetch_orthologs(group_id)
