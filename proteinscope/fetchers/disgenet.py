"""DisGeNET gene-disease association fetcher.

Free tier REST API. Rate limit: conservative 1 req/sec via asyncio.Semaphore + sleep.
No aiolimiter dependency needed.

# Citation: 김박사 (Codex GPT-5.4), Grand Consortium V3 2026-03-20
# API docs: https://disgenet.com/api/
"""

from __future__ import annotations

import asyncio

import httpx

_BASE = "https://www.disgenet.org/api"

# Global semaphore: 1 concurrent request. Combined with asyncio.sleep(1.0)
# after each request this enforces ≤1 req/sec as recommended by 김박사.
_DISGENET_SEM = asyncio.Semaphore(1)


async def fetch_disgenet(gene_symbol: str) -> list[dict]:
    """Return gene-disease associations from DisGeNET for a gene symbol.

    Args:
        gene_symbol: HGNC gene symbol, e.g. "EGFR", "BRCA1".

    Returns:
        List of dicts: [{disease_name, score, disease_id, evidence_count, source}]
        Empty list on any error.
    """
    if not gene_symbol:
        return []

    async with _DISGENET_SEM:
        try:
            url = f"{_BASE}/gda/gene/{gene_symbol}"
            async with httpx.AsyncClient(timeout=30) as client:
                resp = await client.get(
                    url,
                    headers={"accept": "application/json"},
                )
            # Rate limit: 1 req/sec (김박사 recommendation)
            await asyncio.sleep(1.0)

            if resp.status_code == 429:
                # Back off and retry once
                await asyncio.sleep(5.0)
                async with httpx.AsyncClient(timeout=30) as client:
                    resp = await client.get(url, headers={"accept": "application/json"})

            if resp.status_code != 200:
                return []

            data = resp.json()
            if not isinstance(data, list):
                return []

            results = []
            for item in data[:50]:  # cap at 50 to avoid very long payloads
                disease_name = (
                    item.get("disease_name")
                    or item.get("diseaseName")
                    or item.get("disease", {}).get("name", "")
                    or ""
                )
                disease_id = (
                    item.get("diseaseId")
                    or item.get("disease_id")
                    or item.get("disease", {}).get("diseaseId", "")
                    or ""
                )
                score = float(item.get("score") or item.get("gdaScore") or 0.0)
                ev_count = item.get("ei") or item.get("evidenceIndex") or None

                if disease_name:
                    results.append({
                        "disease_name": disease_name,
                        "disease_id": disease_id,
                        "score": score,
                        "evidence_count": int(ev_count) if ev_count is not None else None,
                        "source": "DisGeNET",
                    })
            return results

        except Exception:
            return []
