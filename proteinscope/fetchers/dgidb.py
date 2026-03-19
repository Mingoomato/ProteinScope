"""DGIdb fetcher — drug-gene interactions via GraphQL.

DGIdb aggregates drug-gene interactions from 30+ sources. No API key required.
Used as a secondary source to cross-validate PharmGKB drug interactions.

GraphQL endpoint: https://dgidb.org/api/graphql
"""

from __future__ import annotations

import httpx

_ENDPOINT = "https://dgidb.org/api/graphql"

_QUERY = """
query($genes: [String!]) {
  genes(names: $genes) {
    nodes {
      name
      interactions {
        drug { name conceptId }
        interactionTypes { type directionality }
        interactionScore
        publications { pmid }
        sources { fullName }
      }
    }
  }
}
"""


async def fetch_dgidb_interactions(gene_symbol: str) -> list[dict]:
    """Return drug-gene interactions for a gene symbol from DGIdb.

    Each dict in the returned list has the keys:
      drug, interactionTypes, interactionScore, publications, sources
    """
    async with httpx.AsyncClient() as client:
        try:
            r = await client.post(
                _ENDPOINT,
                json={"query": _QUERY, "variables": {"genes": [gene_symbol]}},
                timeout=30,
            )
            if r.status_code != 200:
                return []
            data = r.json().get("data", {}) or {}
            nodes = (data.get("genes") or {}).get("nodes") or []
            if not nodes:
                return []
            return nodes[0].get("interactions") or []
        except Exception:
            return []
