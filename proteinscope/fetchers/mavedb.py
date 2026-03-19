"""MaveDB deep mutational scanning score-set fetcher.

REST: GET https://api.mavedb.org/api/v1/target-genes/{gene}/score-sets
Returns list of score-set metadata dicts with keys: urn, targetGene, title, numVariants.

Usage::

    from fetchers.mavedb import fetch_mavedb_scores

    score_sets = await fetch_mavedb_scores("BRCA1")
    for ss in score_sets:
        print(ss["urn"], ss["numVariants"])
"""

from __future__ import annotations

import httpx


async def fetch_mavedb_scores(gene_name: str) -> list[dict]:
    """Fetch DMS score-set metadata for a gene from MaveDB.

    GET https://api.mavedb.org/api/v1/target-genes/{gene_name}/score-sets

    Args:
        gene_name: HGNC gene symbol (e.g. "BRCA1", "TP53").

    Returns:
        List of dicts with keys: urn, targetGene, title, numVariants.
        Returns [] on any failure (network error, 404, parse error, etc.).
    """
    try:
        async with httpx.AsyncClient(timeout=20.0) as client:
            r = await client.get(
                f"https://api.mavedb.org/api/v1/target-genes/{gene_name}/score-sets"
            )
            r.raise_for_status()
            raw = r.json()

            # MaveDB may return a list directly or wrap in a results key
            if isinstance(raw, list):
                items = raw
            elif isinstance(raw, dict):
                items = raw.get("results", raw.get("data", []))
                if not isinstance(items, list):
                    items = []
            else:
                return []

            result: list[dict] = []
            for item in items:
                if not isinstance(item, dict):
                    continue
                result.append(
                    {
                        "urn": str(item.get("urn", "") or ""),
                        "targetGene": str(item.get("targetGene", gene_name) or gene_name),
                        "title": str(item.get("title", "") or ""),
                        "numVariants": int(item.get("numVariants", 0) or 0),
                    }
                )
            return result

    except Exception:
        return []
