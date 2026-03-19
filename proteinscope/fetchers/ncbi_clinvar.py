"""NCBI ClinVar fetcher.

Searches ClinVar by gene symbol and extracts clinical significance
and associated disease names from VCV XML records.
"""

from __future__ import annotations

import os
from typing import Optional

import httpx
import xmltodict

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _api_key_params() -> dict:
    key = os.getenv("NCBI_API_KEY")
    return {"api_key": key} if key else {}


async def search_clinvar(gene_symbol: str, client: httpx.AsyncClient) -> list[str]:
    """Search ClinVar for a gene symbol, return list of ClinVar UIDs."""
    params = {
        "db": "clinvar",
        "term": f"{gene_symbol}[gene]",
        "retmax": 20,
        "retmode": "json",
        **_api_key_params(),
    }
    r = await client.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    return data.get("esearchresult", {}).get("idlist", [])


async def fetch_clinvar_vcv(uid: str, client: httpx.AsyncClient) -> dict:
    """Fetch a ClinVar VCV record as parsed XML dict."""
    params = {
        "db": "clinvar",
        "id": uid,
        "rettype": "vcv",
        "retmode": "xml",
        "is_variationid": "true",
        **_api_key_params(),
    }
    r = await client.get(f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    return xmltodict.parse(r.text)


def parse_clinical_entries(vcv_data: dict) -> list[dict]:
    """Extract disease name + significance from a VCV XML parsed dict."""
    entries = []
    try:
        vcv = vcv_data.get("ClinVarResult-Set", {}).get("VariationArchive", {})
        interpretations = vcv.get("InterpretedRecord", {})
        rcvs = interpretations.get("RCVList", {}).get("RCVAccession", [])
        if isinstance(rcvs, dict):
            rcvs = [rcvs]
        for rcv in rcvs:
            classification = rcv.get("Interpretation", {})
            significance = classification.get("Description", "unknown")
            condition = rcv.get("ClassifiedConditionList", {}).get("ClassifiedCondition", {})
            if isinstance(condition, list):
                condition = condition[0]
            disease = condition.get("#text", "") if isinstance(condition, dict) else str(condition)
            if disease:
                entries.append({"disease": disease, "significance": significance, "variant": None})
    except (KeyError, TypeError, AttributeError):
        pass
    return entries


async def fetch_all_clinvar_entries(gene_symbol: str) -> list[dict]:
    """Top-level: search ClinVar for a gene and return all clinical entries."""
    async with httpx.AsyncClient() as client:
        uids = await search_clinvar(gene_symbol, client)
        entries = []
        for uid in uids[:10]:  # cap at 10 to avoid rate limits
            try:
                vcv_data = await fetch_clinvar_vcv(uid, client)
                entries.extend(parse_clinical_entries(vcv_data))
            except Exception:
                continue
        return entries
