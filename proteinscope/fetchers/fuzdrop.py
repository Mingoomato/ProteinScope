"""FuzDrop LLPS propensity fetcher.

REST API: https://fuzdrop.bio.unipd.it/api/predict
No API key needed. Returns LLPS propensity scores per residue.
"""
from __future__ import annotations
import logging
import httpx

_log = logging.getLogger(__name__)
_BASE = "https://fuzdrop.bio.unipd.it/api/predict"

async def fetch_fuzdrop(sequence: str) -> dict:
    """Fetch FuzDrop LLPS phase-separation propensity.

    Returns:
        {
            "llps_score": float,            # overall propensity 0-1
            "per_residue_scores": [float],  # per-residue FuzDrop scores
            "droplet_regions": [{"start": int, "end": int}]
        }
    Returns empty dict on failure.
    """
    if not sequence or len(sequence) < 10:
        return {}
    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            resp = await client.post(
                _BASE,
                json={"sequence": sequence},
                headers={"Content-Type": "application/json"},
            )
            resp.raise_for_status()
            data = resp.json()

            per_residue = data.get("per_residue_scores") or data.get("scores") or []
            llps_score = data.get("llps_score") or data.get("score") or (
                sum(per_residue) / len(per_residue) if per_residue else 0.0
            )

            droplet_regions = data.get("droplet_regions") or data.get("regions") or []
            if not droplet_regions and per_residue:
                # Auto-detect regions with score > 0.6
                in_reg = False
                reg_start = 0
                for i, sc in enumerate(per_residue):
                    if float(sc) > 0.6 and not in_reg:
                        in_reg = True; reg_start = i
                    elif float(sc) <= 0.6 and in_reg:
                        if (i - reg_start) >= 10:
                            droplet_regions.append({"start": reg_start + 1, "end": i})
                        in_reg = False
                if in_reg and (len(per_residue) - reg_start) >= 10:
                    droplet_regions.append({"start": reg_start + 1, "end": len(per_residue)})

            return {
                "llps_score": float(llps_score),
                "per_residue_scores": [float(s) for s in per_residue],
                "droplet_regions": droplet_regions,
            }
    except Exception as exc:
        _log.debug("FuzDrop fetch failed: %s", exc)
        return {}
