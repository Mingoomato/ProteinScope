"""Human Protein Atlas (HPA) REST API fetcher.

Source: Uhlén et al. (2015) Science 347:1260419 — Human Protein Atlas.
API documentation: https://www.proteinatlas.org/about/download
No API key required. Returns empty list on any failure.
"""

from __future__ import annotations

import logging
from typing import Optional

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# HPA field name normalisation map
# ---------------------------------------------------------------------------
# HPA API returns abbreviated column codes; we expand them to human-readable
# field names consistent with the HPA download format described at:
#   https://www.proteinatlas.org/about/download (column definitions)
_TUMOR_FIELD_MAP: dict[str, str] = {
    "g":                 "gene",
    "gs":                "gene_symbol",
    "scl":               "subcellular_location",
    "rnaCancerCategory": "cancer_category",
    "rnaCancerMax":      "tumor_rna_max",
    "rnaNormalMax":      "normal_rna_max",
}

_PROTEIN_FIELD_MAP: dict[str, str] = {
    "g":               "gene",
    "gs":              "gene_symbol",
    "scl":             "subcellular_location",
    "rna_tissue_type": "rna_tissue_type",
    "uniprot":         "uniprot_id",
}


def _normalise(record: dict, field_map: dict[str, str]) -> dict:
    """Return a new dict with HPA column codes replaced by long-form field names."""
    out: dict = {}
    for raw_key, value in record.items():
        mapped = field_map.get(raw_key, raw_key)
        out[mapped] = value
    return out


# ---------------------------------------------------------------------------
# Tumor expression endpoint
# ---------------------------------------------------------------------------

async def fetch_hpa_tumor_expression(tumor_type: str) -> list[dict]:
    """Fetch cancer-vs-normal RNA expression data for a tumor type from HPA.

    Source: Uhlén et al. (2015) Science 347:1260419 — Human Protein Atlas
            tissue and cancer expression atlas.
    API: https://www.proteinatlas.org/about/download

    Args:
        tumor_type: Cancer keyword accepted by HPA search (e.g. "melanoma",
                    "breast cancer", "glioblastoma").

    Returns:
        List of dicts with keys:
            gene, gene_symbol, subcellular_location,
            tumor_rna_max, normal_rna_max, cancer_category
        Returns [] on any network error or malformed response.
    """
    import httpx

    url = (
        "https://www.proteinatlas.org/api/search_download.php"
        f"?search={tumor_type}"
        "&format=json"
        "&columns=g,gs,scl,rnaCancerCategory,rnaCancerMax,rnaNormalMax"
        "&compress=no"
    )

    try:
        async with httpx.AsyncClient(timeout=20.0) as client:
            response = await client.get(url)
            response.raise_for_status()
            raw: list[dict] = response.json()
            if not isinstance(raw, list):
                log.warning("HPA tumor endpoint returned non-list for %s", tumor_type)
                return []
            return [_normalise(rec, _TUMOR_FIELD_MAP) for rec in raw]
    except Exception as exc:
        log.warning("fetch_hpa_tumor_expression failed for %r: %s", tumor_type, exc)
        return []


# ---------------------------------------------------------------------------
# Protein / gene expression endpoint
# ---------------------------------------------------------------------------

async def fetch_hpa_protein_expression(gene_symbol: str) -> dict:
    """Fetch tissue expression profile for a single gene from HPA.

    Source: Uhlén et al. (2015) Science 347:1260419 — Human Protein Atlas
            tissue expression resource.
    API: https://www.proteinatlas.org/about/download

    Args:
        gene_symbol: HGNC gene symbol (e.g. "EGFR", "CD19").

    Returns:
        Dict with keys: gene_symbol, subcellular_location, rna_tissue_type,
        uniprot_id.  Returns {} on any error or if gene is not found.
    """
    import httpx

    url = (
        "https://www.proteinatlas.org/api/search_download.php"
        f"?search={gene_symbol}"
        "&format=json"
        "&columns=g,gs,scl,rna_tissue_type,uniprot"
        "&compress=no"
    )

    try:
        async with httpx.AsyncClient(timeout=20.0) as client:
            response = await client.get(url)
            response.raise_for_status()
            raw = response.json()
            if not isinstance(raw, list) or not raw:
                return {}
            # Return first matching record (most relevant hit for exact gene symbol)
            return _normalise(raw[0], _PROTEIN_FIELD_MAP)
    except Exception as exc:
        log.warning("fetch_hpa_protein_expression failed for %r: %s", gene_symbol, exc)
        return {}
