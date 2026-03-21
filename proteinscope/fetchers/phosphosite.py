"""PhosphoSitePlus PTM site fetcher.

Primary: PhosphoSitePlus public dataset API (reads PHOSPHOSITE_API_KEY env var).
Fallback: Parse UniProt entry["features"] where featureType == "Modified residue".

Returns list of PTM site dicts. Always returns empty list on failure.
"""
from __future__ import annotations
import os
import logging
from typing import Optional
import httpx

_log = logging.getLogger(__name__)
_PSP_BASE = "https://www.phosphosite.org/api/query"

async def fetch_ptm_sites(gene_name: str, uniprot_entry: Optional[dict] = None) -> list[dict]:
    """Fetch PTM sites for a gene from PhosphoSitePlus (+ UniProt fallback).

    Returns list of dicts with keys:
        site (str): e.g. "S473"
        residue (str): single-letter AA code
        position (int): 1-based
        modification_type (str): e.g. "Phosphoserine", "Ubiquitination"
        known_kinases (list[str]): e.g. ["AKT1", "PDK1"]
        functional_effect (str): e.g. "activation", "inhibition", "unknown"
        source (str): "PhosphoSitePlus" or "UniProt"
        pmids (list[str])
    """
    results = []

    # Try PhosphoSitePlus API first
    api_key = os.getenv("PHOSPHOSITE_API_KEY", "")
    if api_key:
        try:
            async with httpx.AsyncClient(timeout=15.0) as client:
                resp = await client.get(
                    f"{_PSP_BASE}/ptm",
                    params={"gene": gene_name, "organism": "human"},
                    headers={"Authorization": f"Bearer {api_key}"},
                )
                if resp.status_code == 200:
                    data = resp.json()
                    for site in (data.get("sites") or []):
                        results.append({
                            "site": site.get("mod_site", ""),
                            "residue": site.get("residue", ""),
                            "position": int(site.get("position", 0)),
                            "modification_type": site.get("modification", ""),
                            "known_kinases": site.get("kinases", []),
                            "functional_effect": site.get("effect", "unknown"),
                            "source": "PhosphoSitePlus",
                            "pmids": site.get("pmids", []),
                        })
        except Exception as exc:
            _log.debug("PhosphoSitePlus API failed for %s: %s", gene_name, exc)

    # Fallback: UniProt entry["features"] where featureType == "Modified residue"
    if not results and uniprot_entry:
        try:
            for feat in (uniprot_entry.get("features") or []):
                if feat.get("type") != "Modified residue":
                    continue
                loc = feat.get("location", {})
                start = loc.get("start", {}).get("value", 0)
                desc = feat.get("description", "")
                # Try to extract AA from sequence at position
                residue = ""
                try:
                    seq = uniprot_entry.get("sequence", {}).get("sequence", "")
                    if seq and 0 < start <= len(seq):
                        residue = seq[start - 1]
                except Exception:
                    pass
                results.append({
                    "site": f"{residue}{start}",
                    "residue": residue,
                    "position": start,
                    "modification_type": desc,
                    "known_kinases": [],
                    "functional_effect": "unknown",
                    "source": "UniProt",
                    "pmids": [],
                })
        except Exception as exc:
            _log.debug("UniProt PTM fallback failed for %s: %s", gene_name, exc)

    return results
