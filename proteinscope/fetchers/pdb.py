"""RCSB PDB fetcher for experimental 3D structures.

Uses the RCSB Search API (POST, JSON) to find structures containing
the given UniProt accession as a polymer entity.
No API key required.
"""

from __future__ import annotations

import httpx

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA = "https://data.rcsb.org/rest/v1/core/entry"


async def search_pdb_by_uniprot(accession: str) -> list[str]:
    """Return list of PDB IDs that contain the given UniProt accession."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession,
            },
        },
        "return_type": "entry",
        "request_options": {"results_verbosity": "minimal", "paginate": {"start": 0, "rows": 20}},
    }
    async with httpx.AsyncClient(timeout=30) as client:
        r = await client.post(RCSB_SEARCH, json=query)
        if r.status_code == 204:
            return []
        r.raise_for_status()
        data = r.json()
        return [hit["identifier"] for hit in data.get("result_set", [])]


async def fetch_pdb_entry(pdb_id: str, client: httpx.AsyncClient) -> dict:
    """Fetch entry metadata (resolution, method) from RCSB Data API."""
    r = await client.get(f"{RCSB_DATA}/{pdb_id}", timeout=30)
    r.raise_for_status()
    return r.json()


async def fetch_experimental_structures(accession: str) -> list[dict]:
    """Return list of dicts: {pdb_id, resolution, method} for a UniProt accession."""
    pdb_ids = await search_pdb_by_uniprot(accession)
    structures = []
    async with httpx.AsyncClient(timeout=30) as client:
        for pdb_id in pdb_ids:
            try:
                entry = await fetch_pdb_entry(pdb_id, client)
                expt = entry.get("exptl", [{}])[0]
                method = expt.get("method", "")
                refine = entry.get("refine", [{}])
                resolution = None
                if refine:
                    resolution = refine[0].get("ls_d_res_high")
                structures.append({
                    "pdb_id": pdb_id,
                    "resolution": resolution,
                    "method": method,
                })
            except Exception:
                structures.append({"pdb_id": pdb_id, "resolution": None, "method": ""})
    return structures
