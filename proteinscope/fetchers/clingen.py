"""ClinGen Gene-Disease Validity fetcher.

Returns Definitive/Strong/Moderate/Limited classifications for gene-disease pairs.
Free REST API, no key required.

# Citation: 노박사 (Gemini 2.5 Pro), Grand Consortium V3 2026-03-20
# API: https://search.clinicalgenome.org/kb/gene-validity
# Docs: https://clinicalgenome.org/curation-activities/gene-disease-validity/
"""

from __future__ import annotations

import httpx

_BASE = "https://search.clinicalgenome.org/kb/gene-validity"

# ClinGen classification strength order (for filtering)
CLASSIFICATION_ORDER = {
    "Definitive": 5,
    "Strong": 4,
    "Moderate": 3,
    "Limited": 2,
    "Disputed": 1,
    "Refuted": 0,
    "No Known Disease Relationship": 0,
}


async def fetch_clingen(gene_symbol: str) -> list[dict]:
    """Return ClinGen gene-disease validity classifications for a gene.

    Args:
        gene_symbol: HGNC gene symbol, e.g. "BRCA1", "CFTR".

    Returns:
        List of dicts: [{disease_name, classification, mondo_id, pmids}]
        Empty list on any error.
    """
    if not gene_symbol:
        return []

    try:
        params = {
            "search": gene_symbol,
            "format": "json",
            "limit": 50,
            "skip": 0,
        }
        async with httpx.AsyncClient(timeout=30) as client:
            resp = await client.get(_BASE, params=params)

        if resp.status_code != 200:
            return []

        data = resp.json()
        # ClinGen returns {"gene_validity_list": [...]} or similar
        items = (
            data.get("gene_validity_list")
            or data.get("data")
            or (data if isinstance(data, list) else [])
        )

        results = []
        for item in items:
            # Gene symbol filter — ClinGen search may return partial matches
            item_gene = (
                item.get("gene", {}).get("gene_symbol")
                or item.get("geneSymbol")
                or item.get("gene_symbol")
                or ""
            )
            if item_gene.upper() != gene_symbol.upper():
                continue

            disease_name = (
                item.get("disease", {}).get("label")
                or item.get("disease_label")
                or item.get("diseaseName")
                or ""
            )
            classification = (
                item.get("classification", {}).get("label")
                or item.get("classification")
                or ""
            )
            mondo_id = (
                item.get("disease", {}).get("curie")
                or item.get("disease_id")
                or ""
            )
            pmids = item.get("pmids") or []

            if disease_name and classification:
                results.append({
                    "disease_name": disease_name,
                    "classification": classification,
                    "mondo_id": mondo_id,
                    "pmids": pmids,
                    "source": "ClinGen",
                    "classification_rank": CLASSIFICATION_ORDER.get(classification, 0),
                })

        # Sort by classification strength descending
        results.sort(key=lambda x: x["classification_rank"], reverse=True)
        return results

    except Exception:
        return []
