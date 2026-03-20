"""Human Phenotype Ontology (HPO) fetcher.

Returns HPO terms associated with a gene via NCBI gene ID.
Free REST API, no key required.

# Citation: Grand Consortium V3 2026-03-20
# API: https://hpo.jax.org/api/hpo/gene/{ncbi_gene_id}
# Docs: https://hpo.jax.org/
"""

from __future__ import annotations

import httpx

_BASE = "https://hpo.jax.org/api/hpo/gene"


async def fetch_hpo_terms(ncbi_gene_id: str) -> list[dict]:
    """Return HPO phenotype terms for an NCBI gene ID.

    Args:
        ncbi_gene_id: NCBI Entrez gene ID as string, e.g. "1956" (EGFR).

    Returns:
        List of dicts: [{hpo_id, hpo_name, definition}]
        Empty list on any error.
    """
    if not ncbi_gene_id:
        return []

    try:
        url = f"{_BASE}/{ncbi_gene_id}"
        async with httpx.AsyncClient(timeout=30) as client:
            resp = await client.get(url, headers={"accept": "application/json"})

        if resp.status_code != 200:
            return []

        data = resp.json()
        # HPO API returns {"termCount": N, "terms": [{ontologyId, name, definition, ...}]}
        terms = data.get("terms") or data.get("associations") or []

        results = []
        for term in terms:
            hpo_id = term.get("ontologyId") or term.get("hpoId") or ""
            hpo_name = term.get("name") or ""
            definition = term.get("definition") or ""
            if hpo_id and hpo_name:
                results.append({
                    "hpo_id": hpo_id,
                    "hpo_name": hpo_name,
                    "definition": definition,
                })
        return results

    except Exception:
        return []
