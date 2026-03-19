"""ADMET property fetcher using the pkCSM REST API.

pkCSM (Pires et al. 2015) provides in-silico ADMET predictions from SMILES.
When the pkCSM server is unavailable, Lipinski rule-of-5 violations are computed
locally from the SMILES string via atom-counting heuristics.

Sources:
  - pkCSM API  https://biosig.lab.uq.edu.au/pkcsm/api/predictions
  - ChEMBL API https://www.ebi.ac.uk/chembl/api/data (SMILES lookup fallback)

Reference:
  # pkCSM: Pires DE 2015 J Med Chem doi:10.1021/acs.jmedchem.4c02234
  # Lipinski rule of 5: Lipinski CA 2001 Adv Drug Deliv Rev
  # doi:10.1016/S0169-409X(00)00129-0
"""

from __future__ import annotations

import asyncio
import logging
import re
from typing import Optional

import httpx

_log = logging.getLogger(__name__)

# Delay between pkCSM requests (seconds) to respect rate limit
_PKCSM_RATE_DELAY = 0.3

_PKCSM_URL   = "https://biosig.lab.uq.edu.au/pkcsm/api/predictions"
_CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


# ---------------------------------------------------------------------------
# Lipinski rule-of-5 (local heuristic)
# ---------------------------------------------------------------------------

async def _compute_lipinski(smiles: str) -> int:
    """Count Lipinski rule-of-5 violations from a SMILES string.

    Uses atom-counting heuristics derived directly from the SMILES notation:
      - MW  proxy: 14*nC + 16*nO + 14*nN + 32*nS + 1*(explicit H implied)
      - HBA: count of N and O atoms
      - HBD: count of NH and OH groups (N or O followed by H in SMILES)
      - LogP proxy: (nC - nO - nN) / max(nC, 1)  [very rough]

    # Lipinski rule of 5: Lipinski CA 2001 Adv Drug Deliv Rev
    # doi:10.1016/S0169-409X(00)00129-0
    # MW > 500, HBD > 5, HBA > 10, LogP > 5 trigger violations.

    Args:
        smiles: SMILES string.

    Returns:
        Integer count of violated rules (0–4).
    """
    if not smiles:
        return 0

    try:
        s = smiles.upper()

        # Count heavy atom types from SMILES (case-insensitive atomic symbols)
        n_c  = len(re.findall(r"C(?!L|O|A)", s))          # C but not Cl, Co, Ca
        n_o  = len(re.findall(r"O", s))
        n_n  = len(re.findall(r"N", s))
        n_s  = len(re.findall(r"S(?!E|I)", s))            # S but not Se, Si
        n_cl = len(re.findall(r"CL", s))
        n_br = len(re.findall(r"BR", s))
        n_f  = len(re.findall(r"F", s))

        # Molecular weight proxy (Da)
        mw_proxy = (
            12 * n_c
            + 16 * n_o
            + 14 * n_n
            + 32 * n_s
            + 35 * n_cl
            + 80 * n_br
            + 19 * n_f
        )

        # HBA: number of N + O atoms
        hba = n_n + n_o

        # HBD: NH or OH groups — look for [NH], [OH], or N/O adjacent to H in SMILES
        # In SMILES bracket notation [NH2], [NH], [OH]; also bare N (amine implicit H)
        hbd = len(re.findall(r"\[(?:NH2?|OH)\]", smiles, re.IGNORECASE))
        # Also count bare nitrogen atoms that are likely amines (not aromatic, no bracket)
        # This is approximate; full HBD requires a proper parser
        hbd += len(re.findall(r"(?<!\[)N(?!H)(?![a-z])", smiles))

        # LogP proxy: Crippen-inspired heuristic (very rough)
        logp_proxy = (n_c - n_o - n_n + 0.5 * n_cl + 1.0 * n_br) / max(n_c, 1) * 2.5

        violations = 0
        if mw_proxy > 500:   violations += 1   # MW > 500
        if hbd      > 5:     violations += 1   # HBD > 5
        if hba      > 10:    violations += 1   # HBA > 10
        if logp_proxy > 5:   violations += 1   # LogP > 5

        return violations
    except Exception as exc:
        _log.debug("_compute_lipinski failed for smiles=%s: %s", smiles[:20], exc)
        return 0


# ---------------------------------------------------------------------------
# ChEMBL SMILES lookup
# ---------------------------------------------------------------------------

async def _fetch_smiles_from_chembl(drug_name: str, client: httpx.AsyncClient) -> str:
    """Attempt to retrieve a canonical SMILES from ChEMBL by drug name.

    Args:
        drug_name: Drug name string (approximate match attempted).
        client:    Shared httpx.AsyncClient instance.

    Returns:
        SMILES string, or empty string if not found.
    """
    try:
        r = await client.get(
            f"{_CHEMBL_BASE}/molecule",
            params={
                "pref_name__iexact": drug_name,
                "format":            "json",
                "limit":             1,
            },
            timeout=15.0,
        )
        r.raise_for_status()
        molecules = r.json().get("molecules", [])
        if not molecules:
            return ""
        struct = molecules[0].get("molecule_structures") or {}
        return str(struct.get("canonical_smiles", ""))
    except Exception as exc:
        _log.debug("ChEMBL SMILES lookup failed for %s: %s", drug_name, exc)
        return ""


# ---------------------------------------------------------------------------
# pkCSM prediction
# ---------------------------------------------------------------------------

async def _call_pkcsm(smiles: str, client: httpx.AsyncClient) -> dict:
    """POST a single SMILES to the pkCSM prediction endpoint.

    # pkCSM: Pires DE 2015 J Med Chem doi:10.1021/acs.jmedchem.4c02234

    Args:
        smiles: Canonical SMILES string.
        client: Shared httpx.AsyncClient instance.

    Returns:
        Raw response JSON dict from pkCSM, or {} on failure.
    """
    try:
        payload = {"smiles": smiles, "input": [smiles]}
        r = await client.post(_PKCSM_URL, json=payload, timeout=30.0)
        r.raise_for_status()
        return r.json()
    except Exception as exc:
        _log.debug("pkCSM call failed for smiles=%s: %s", smiles[:20], exc)
        return {}


def _extract_pkcsm_fields(raw: dict) -> dict:
    """Parse pkCSM API response into a flat dict of ADMET properties.

    pkCSM returns a JSON structure; field names vary by API version.
    We attempt multiple plausible key names and fall back to None.

    Args:
        raw: Raw pkCSM response dict.

    Returns:
        Dict with keys: hia, bbb_permeant, cyp3a4_inhibitor, herg_blocker,
        ames_toxic.
    """
    # pkCSM v2 response may nest results under "predictions" or be flat
    data = raw.get("predictions", raw) if isinstance(raw, dict) else {}
    if isinstance(data, list) and data:
        data = data[0]
    if not isinstance(data, dict):
        return {}

    def _bool(key: str) -> Optional[bool]:
        val = data.get(key)
        if val is None:
            return None
        if isinstance(val, bool):
            return val
        if isinstance(val, (int, float)):
            return bool(val)
        # Some endpoints return "YES"/"NO" or "Positive"/"Negative"
        return str(val).strip().lower() in ("yes", "positive", "true", "1")

    def _float(key: str) -> Optional[float]:
        val = data.get(key)
        if val is None:
            return None
        try:
            return float(val)
        except (TypeError, ValueError):
            return None

    from typing import Optional  # local import to satisfy type checker

    return {
        # Human intestinal absorption — pkCSM calls this "Absorption_HIA_Human"
        "hia":             _float("Absorption_HIA_Human") or _float("hia"),
        # BBB permeability
        "bbb_permeant":    _bool("Distribution_BBB") or _bool("bbb_permeant"),
        # CYP3A4 inhibition
        "cyp3a4_inhibitor": _bool("Metabolism_CYP3A4_inhibitor") or _bool("cyp3a4_inhibitor"),
        # hERG cardiotoxicity
        # hERG cardiac risk: Sanguinetti MC 2006 Nature doi:10.1038/nature04710
        "herg_blocker":    _bool("Toxicity_hERG_I") or _bool("Toxicity_hERG_II") or _bool("herg_blocker"),
        # Ames mutagenicity
        "ames_toxic":      _bool("Toxicity_AMES") or _bool("ames_toxic"),
    }


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

async def fetch_admet_profiles(drug_list: list[dict]) -> list[dict]:
    """Fetch ADMET profiles for a list of drug dicts via pkCSM.

    For each drug that lacks a SMILES key, attempts to retrieve it from
    ChEMBL by name.  Profiles are fetched sequentially with a 0.3 s inter-
    request delay to respect the pkCSM rate limit.  Partial results are
    returned if the server fails for a subset of drugs.

    # pkCSM: Pires DE 2015 J Med Chem doi:10.1021/acs.jmedchem.4c02234
    # Lipinski rule of 5: Lipinski CA 2001 Adv Drug Deliv Rev
    # doi:10.1016/S0169-409X(00)00129-0

    Args:
        drug_list: List of dicts, each expected to contain at least one of
                   'smiles' or 'drug_name'.  'drug_id' is used as an opaque
                   identifier; if absent a sequential index is substituted.

    Returns:
        List of ADMET profile dicts with keys:
            drug_id, drug_name, smiles, hia, bbb_permeant, cyp3a4_inhibitor,
            herg_blocker, ames_toxic, lipinski_violations.
        Empty list if all requests fail or drug_list is empty.
    """
    if not drug_list:
        return []

    results: list[dict] = []

    async with httpx.AsyncClient(timeout=30.0) as client:
        for idx, drug in enumerate(drug_list):
            try:
                drug_id   = str(drug.get("drug_id")   or drug.get("id")   or f"drug_{idx}")
                drug_name = str(drug.get("drug_name") or drug.get("name") or drug_id)
                smiles    = str(drug.get("smiles", "") or "").strip()

                # Attempt ChEMBL SMILES lookup if SMILES is missing
                if not smiles and drug_name:
                    smiles = await _fetch_smiles_from_chembl(drug_name, client)

                # Compute local Lipinski violations (always available)
                lipinski_violations = await _compute_lipinski(smiles)

                # Call pkCSM if SMILES is available
                admet_fields: dict = {
                    "hia":             None,
                    "bbb_permeant":    None,
                    "cyp3a4_inhibitor": None,
                    "herg_blocker":    None,
                    "ames_toxic":      None,
                }
                if smiles:
                    raw_response = await _call_pkcsm(smiles, client)
                    if raw_response:
                        admet_fields = _extract_pkcsm_fields(raw_response)

                results.append({
                    "drug_id":             drug_id,
                    "drug_name":           drug_name,
                    "smiles":              smiles or None,
                    "hia":                 admet_fields.get("hia"),
                    "bbb_permeant":        admet_fields.get("bbb_permeant"),
                    "cyp3a4_inhibitor":    admet_fields.get("cyp3a4_inhibitor"),
                    "herg_blocker":        admet_fields.get("herg_blocker"),
                    "ames_toxic":          admet_fields.get("ames_toxic"),
                    "lipinski_violations": lipinski_violations,
                })

                # Rate-limit between requests
                if idx < len(drug_list) - 1:
                    await asyncio.sleep(_PKCSM_RATE_DELAY)

            except Exception as exc:
                _log.debug("fetch_admet_profiles failed for drug[%d]: %s", idx, exc)
                # Append a minimal record so the downstream analyzer can still work
                results.append({
                    "drug_id":             str(drug.get("drug_id") or f"drug_{idx}"),
                    "drug_name":           str(drug.get("drug_name") or f"drug_{idx}"),
                    "smiles":              None,
                    "hia":                 None,
                    "bbb_permeant":        None,
                    "cyp3a4_inhibitor":    None,
                    "herg_blocker":        None,
                    "ames_toxic":          None,
                    "lipinski_violations": 0,
                })

    return results
