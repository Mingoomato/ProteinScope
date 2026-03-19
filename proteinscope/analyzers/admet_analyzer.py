"""ADMET Analyzer — drug ADMET risk profiling and clinical liability synthesis.

Merges DrugInteraction objects with pkCSM-derived ADMET properties, assigns a
three-tier overall risk classification (low / medium / high), and synthesises a
Gemini-powered clinical liability summary for the set of drugs interacting with
the query protein.

Key citations embedded at threshold/formula sites:
  # hERG cardiac risk: Sanguinetti MC 2006 Nature doi:10.1038/nature04710
  # Ames mutagenicity: Maron DM 1983 Mutat Res doi:10.1016/0027-5107(83)90283-0
  # Lipinski rule of 5: Lipinski CA 2001 Adv Drug Deliv Rev
  # doi:10.1016/S0169-409X(00)00129-0
"""

from __future__ import annotations

import json as _json
import logging
from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel, Field

from core.evidence import DataProvenance, EvidenceGrade

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class DrugADMET(BaseModel):
    """ADMET profile for a single drug."""

    drug_name:          str
    drug_id:            str
    smiles:             Optional[str]   = None
    hia:                Optional[float] = None       # human intestinal absorption 0–1
    bbb_permeant:       Optional[bool]  = None
    cyp3a4_inhibitor:   Optional[bool]  = None
    herg_blocker:       Optional[bool]  = None
    ames_toxic:         Optional[bool]  = None
    lipinski_violations: int            = 0
    # "low" | "medium" | "high"
    overall_risk:       str             = "unknown"
    provenance:         Optional[DataProvenance] = None


class ADMETReport(BaseModel):
    """Full ADMET analysis report for a gene and its associated drugs."""

    gene:            str
    drug_profiles:   List[DrugADMET]
    high_risk_count: int
    key_liabilities: List[str]
    gemini_summary:  str = ""
    timestamp:       datetime


# ---------------------------------------------------------------------------
# Risk classification
# ---------------------------------------------------------------------------

def _classify_risk(
    herg_blocker:    Optional[bool],
    ames_toxic:      Optional[bool],
    lipinski_violations: int,
    cyp3a4_inhibitor: Optional[bool],
) -> str:
    """Assign a three-tier overall ADMET risk category.

    High risk if:
      - hERG blockade (cardiac QT prolongation risk)
        # hERG cardiac risk: Sanguinetti MC 2006 Nature doi:10.1038/nature04710
      - Ames mutagenicity positive
        # Ames test: Maron DM 1983 Mutat Res doi:10.1016/0027-5107(83)90283-0

    Medium risk if:
      - >= 2 Lipinski violations (drug-likeness concern)
        # Lipinski rule of 5: Lipinski CA 2001 Adv Drug Deliv Rev
        # doi:10.1016/S0169-409X(00)00129-0
      - CYP3A4 inhibition (DDI potential)

    Otherwise: low risk.

    Args:
        herg_blocker:        hERG cardiac liability flag.
        ames_toxic:          Ames mutagenicity flag.
        lipinski_violations: Number of Lipinski rule-of-5 violations.
        cyp3a4_inhibitor:    CYP3A4 metabolic inhibition flag.

    Returns:
        "high", "medium", or "low".
    """
    if herg_blocker or ames_toxic:
        return "high"
    if lipinski_violations >= 2 or cyp3a4_inhibitor:
        return "medium"
    return "low"


# ---------------------------------------------------------------------------
# Gemini synthesis
# ---------------------------------------------------------------------------

async def _synthesize_admet(
    gene: str,
    profiles: list[DrugADMET],
    high_risk_count: int,
) -> str:
    """Call Gemini for a clinical ADMET liability narrative.

    Args:
        gene:            Gene symbol.
        profiles:        List of DrugADMET objects.
        high_risk_count: Pre-computed count of high-risk drugs.

    Returns:
        Gemini narrative string; empty string on failure.
    """
    try:
        from core.gemini_interpreter import _call

        profiles_json = _json.dumps(
            [
                {
                    "drug_name":            p.drug_name,
                    "hia":                  p.hia,
                    "bbb_permeant":         p.bbb_permeant,
                    "cyp3a4_inhibitor":     p.cyp3a4_inhibitor,
                    "herg_blocker":         p.herg_blocker,
                    "ames_toxic":           p.ames_toxic,
                    "lipinski_violations":  p.lipinski_violations,
                    "overall_risk":         p.overall_risk,
                }
                for p in profiles[:25]
            ],
            indent=2,
        )

        prompt = (
            f"You are a clinical pharmacology and drug-safety expert. "
            f"The following ADMET profiles are for drugs interacting with {gene}:\n"
            f"{profiles_json}\n\n"
            f"High-risk drug count: {high_risk_count} of {len(profiles)}.\n\n"
            "Provide a concise clinical liability summary covering:\n"
            "1. The most clinically significant ADMET liabilities across this drug set.\n"
            "2. hERG-related cardiac risk (QT prolongation) if herg_blocker is present.\n"
            "   # hERG cardiac risk: Sanguinetti MC 2006 Nature doi:10.1038/nature04710\n"
            "3. CYP3A4-mediated drug-drug interaction risk.\n"
            "4. Genotoxicity concerns from Ames-positive compounds.\n"
            "5. Recommendations for medicinal chemistry optimization.\n\n"
            "Return ONLY a plain text paragraph (no JSON, no markdown headers)."
        )

        raw = await _call(prompt)
        return raw.strip() if raw else ""
    except Exception as exc:
        _log.debug("Gemini ADMET synthesis failed for %s: %s", gene, exc)
        return ""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _extract_drug_id(obj) -> str:
    """Extract a drug identifier from a DrugInteraction object or dict."""
    if isinstance(obj, dict):
        return str(
            obj.get("drug_id")
            or obj.get("chembl_id")
            or obj.get("id")
            or ""
        ).strip()
    # Pydantic model or dataclass attribute access
    for attr in ("drug_id", "chembl_id", "id"):
        val = getattr(obj, attr, None)
        if val:
            return str(val).strip()
    return ""


def _extract_drug_name(obj) -> str:
    """Extract a human-readable drug name from a DrugInteraction object or dict."""
    if isinstance(obj, dict):
        return str(
            obj.get("drug_name")
            or obj.get("name")
            or obj.get("drug_id")
            or ""
        ).strip()
    for attr in ("drug_name", "name", "drug_id"):
        val = getattr(obj, attr, None)
        if val:
            return str(val).strip()
    return ""


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_admet_analysis(
    gene: str,
    drug_interactions: list,   # DrugInteraction objects or dicts
    admet_data: list,          # from fetch_admet_profiles()
    step_cb=None,
) -> Optional[ADMETReport]:
    """Run ADMET analysis for drugs associated with a gene.

    Args:
        gene:              Gene symbol, e.g. "CYP3A4".
        drug_interactions: List of DrugInteraction objects or plain dicts
                           describing drugs interacting with the gene.
        admet_data:        Output of fetchers.admet.fetch_admet_profiles().
        step_cb:           Optional async progress callback (receives str).

    Returns:
        ADMETReport, or None if a fatal error occurs.
    """
    async def _step(msg: str) -> None:
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    try:
        # ── 1. Parse drug_interactions into (drug_id → name) map ─────────────
        await _step("[1/4] Parsing drug interactions...")
        interaction_map: dict[str, str] = {}   # drug_id → drug_name
        for di in drug_interactions or []:
            did  = _extract_drug_id(di)
            name = _extract_drug_name(di)
            if did:
                interaction_map[did] = name or did

        # ── 2. Merge with admet_data by drug_id ───────────────────────────────
        await _step("[2/4] Merging ADMET data with drug interaction profiles...")
        admet_by_id: dict[str, dict] = {}
        for rec in admet_data or []:
            did = str(rec.get("drug_id", "")).strip()
            if did:
                admet_by_id[did] = rec

        # Build the union of all known drug IDs
        all_ids = set(interaction_map.keys()) | set(admet_by_id.keys())

        profiles: list[DrugADMET] = []
        for did in all_ids:
            rec       = admet_by_id.get(did, {})
            drug_name = interaction_map.get(did) or rec.get("drug_name") or did

            herg     = rec.get("herg_blocker")
            ames     = rec.get("ames_toxic")
            lip_viol = int(rec.get("lipinski_violations", 0))
            cyp3a4   = rec.get("cyp3a4_inhibitor")

            # ── 3. Compute overall_risk ───────────────────────────────────────
            risk = _classify_risk(herg, ames, lip_viol, cyp3a4)

            provenance = DataProvenance(
                source="pkCSM" if rec else "DrugInteraction list",
                evidence_grade=(
                    EvidenceGrade.COMPUTATIONAL if rec
                    else EvidenceGrade.LITERATURE
                ),
                scientific_caveat=(
                    "pkCSM in-silico ADMET prediction; experimental validation required."
                    if rec else
                    "Drug identity from interaction database; ADMET data unavailable."
                ),
                method="pkCSM REST API" if rec else None,
            )

            profiles.append(DrugADMET(
                drug_name=drug_name,
                drug_id=did,
                smiles=rec.get("smiles"),
                hia=rec.get("hia"),
                bbb_permeant=rec.get("bbb_permeant"),
                cyp3a4_inhibitor=cyp3a4,
                herg_blocker=herg,
                ames_toxic=ames,
                lipinski_violations=lip_viol,
                overall_risk=risk,
                provenance=provenance,
            ))

        # Sort: high risk first, then medium, then low/unknown
        _risk_order = {"high": 0, "medium": 1, "low": 2, "unknown": 3}
        profiles.sort(key=lambda p: _risk_order.get(p.overall_risk, 3))

        high_risk_count = sum(1 for p in profiles if p.overall_risk == "high")

        # Derive key liabilities list
        key_liabilities: list[str] = []
        if any(p.herg_blocker for p in profiles):
            # hERG cardiac risk: Sanguinetti MC 2006 Nature doi:10.1038/nature04710
            herg_names = [p.drug_name for p in profiles if p.herg_blocker]
            key_liabilities.append(
                f"hERG cardiac risk (QT prolongation): {', '.join(herg_names[:5])}"
            )
        if any(p.ames_toxic for p in profiles):
            ames_names = [p.drug_name for p in profiles if p.ames_toxic]
            key_liabilities.append(
                f"Ames mutagenicity: {', '.join(ames_names[:5])}"
            )
        if any(p.cyp3a4_inhibitor for p in profiles):
            cyp_names = [p.drug_name for p in profiles if p.cyp3a4_inhibitor]
            key_liabilities.append(
                f"CYP3A4 inhibition / DDI risk: {', '.join(cyp_names[:5])}"
            )
        high_lip = [p.drug_name for p in profiles if p.lipinski_violations >= 2]
        if high_lip:
            # Lipinski rule of 5: Lipinski CA 2001 Adv Drug Deliv Rev
            # doi:10.1016/S0169-409X(00)00129-0
            key_liabilities.append(
                f"Poor drug-likeness (>=2 Lipinski violations): {', '.join(high_lip[:5])}"
            )

        # ── 4. Gemini synthesis ───────────────────────────────────────────────
        await _step("[4/4] Gemini clinical liability synthesis...")
        gemini_summary = await _synthesize_admet(gene, profiles, high_risk_count)

        return ADMETReport(
            gene=gene,
            drug_profiles=profiles,
            high_risk_count=high_risk_count,
            key_liabilities=key_liabilities,
            gemini_summary=gemini_summary,
            timestamp=datetime.utcnow(),
        )

    except Exception as exc:
        _log.debug("run_admet_analysis failed for %s: %s", gene, exc)
        return None
