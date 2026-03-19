"""PhaSepDB protein phase-separation database fetcher.

PhaSepDB is a manually curated database of phase-separating proteins and their
condensate properties. It provides LLPS propensity scores, condensate type
classifications, and known interaction partners for proteins implicated in
liquid-liquid phase separation (LLPS) and biomolecular condensate formation.

REST endpoint:
    GET https://db.phasep.pro/api/protein?uniprot={acc}

Returns a dict with keys:
    llps_propensity  : float  — numeric LLPS propensity score (0–1)
    condensate_types : list   — e.g. ["stress granule", "P-body"]
    partner_proteins : list   — UniProt accessions of known partners

Falls back to an empty dict {} on any network or parse error.

Reference:
    You S et al. (2020) PhaSepDB: a database of liquid-liquid phase separation
    related proteins. Nucleic Acids Res 48(D1):D354-D359.
    doi:10.1093/nar/gkz927
"""

from __future__ import annotations

from typing import Optional

import httpx

_PHASEPDB_BASE = "https://db.phasep.pro/api/protein"


async def fetch_phasepdb(uniprot_acc: str) -> dict:
    """Fetch phase-separation data for a UniProt accession from PhaSepDB.

    Args:
        uniprot_acc: UniProt accession string, e.g. "P04637" (TP53).

    Returns:
        A dict with keys ``llps_propensity`` (float), ``condensate_types``
        (list[str]), and ``partner_proteins`` (list[str]).  Returns an empty
        dict on any failure so callers never need to guard against None.
    """
    if not uniprot_acc or not uniprot_acc.strip():
        return {}
    try:
        url = f"{_PHASEPDB_BASE}?uniprot={uniprot_acc.strip()}"
        async with httpx.AsyncClient(timeout=20.0) as client:
            r = await client.get(url)
            r.raise_for_status()
            data = r.json()
        # Normalise: guarantee the expected keys are always present
        return {
            "llps_propensity": float(data.get("llps_propensity") or 0.0),
            "condensate_types": list(data.get("condensate_types") or []),
            "partner_proteins": list(data.get("partner_proteins") or []),
        }
    except Exception:
        return {}
