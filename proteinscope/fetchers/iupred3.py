"""IUPred3 disorder prediction fetcher.

REST API: https://iupred3.elte.hu/iupred3API
No API key needed. Returns per-residue disorder scores.
"""
from __future__ import annotations
import logging
import httpx

_log = logging.getLogger(__name__)
_BASE = "https://iupred3.elte.hu/iupred3API"

async def fetch_iupred_scores(sequence: str) -> dict:
    """Fetch per-residue IUPred3 disorder scores.

    Returns:
        {
            "scores": [float per residue, 0.0=ordered, 1.0=disordered],
            "regions": [{"start": int, "end": int, "type": "disordered"}]
        }
    Returns empty dict on failure.
    """
    if not sequence or len(sequence) < 5:
        return {}
    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            resp = await client.get(
                _BASE,
                params={"sequence": sequence, "type": "long"},
            )
            resp.raise_for_status()
            data = resp.json()

            # IUPred3 returns {"iupred2": [...scores...]} or similar
            scores = data.get("iupred2") or data.get("iupred3") or data.get("scores") or []
            if not scores and isinstance(data, list):
                scores = data

            scores = [float(s) for s in scores]

            # Identify disordered regions (score > 0.5 for ≥20 consecutive residues)
            regions = []
            in_region = False
            start = 0
            for i, sc in enumerate(scores):
                if sc > 0.5 and not in_region:
                    in_region = True
                    start = i
                elif sc <= 0.5 and in_region:
                    if (i - start) >= 20:
                        regions.append({"start": start + 1, "end": i, "type": "disordered"})
                    in_region = False
            if in_region and (len(scores) - start) >= 20:
                regions.append({"start": start + 1, "end": len(scores), "type": "disordered"})

            return {"scores": scores, "regions": regions}
    except Exception as exc:
        _log.debug("IUPred3 fetch failed: %s", exc)
        return {}
