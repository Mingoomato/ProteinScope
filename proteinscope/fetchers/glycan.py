"""Glycan analysis data fetcher.

Combines GlyConnect (known glycan compositions from experiment) and a
NetNGlyc-style rule-based fallback for N-glycosylation site prediction
when the DTU server is unreachable.

Sources:
  - GlyConnect REST API  (https://glyconnect.expasy.org)
  - GlyTouCan REST API   (https://api.glytoucan.org)
  - DTU NetNGlyc server  (https://services.healthtech.dtu.dk)
"""

from __future__ import annotations

import asyncio
import logging
import re

import httpx

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# GlyConnect
# ---------------------------------------------------------------------------

async def fetch_glyconnect(uniprot_acc: str) -> list[dict]:
    """Return known glycan compositions for a UniProt accession from GlyConnect.

    Args:
        uniprot_acc: UniProt accession, e.g. "P12345".

    Returns:
        List of dicts with keys: position, glycan_composition, evidence_type,
        glytoucan_id.  Empty list on any error.
    """
    try:
        url = f"https://glyconnect.expasy.org/api/proteins/{uniprot_acc}/compositions"
        async with httpx.AsyncClient(timeout=20.0) as client:
            r = await client.get(url)
            r.raise_for_status()
            raw = r.json()

        results: list[dict] = []
        entries = raw if isinstance(raw, list) else raw.get("compositions", [])
        for entry in entries:
            results.append({
                "position":          entry.get("position"),
                "glycan_composition": entry.get("composition") or entry.get("glycan_composition", ""),
                "evidence_type":     entry.get("evidence_type", ""),
                "glytoucan_id":      entry.get("glytoucan_id") or entry.get("accession", ""),
            })
        return results
    except Exception as exc:
        _log.debug("fetch_glyconnect failed for %s: %s", uniprot_acc, exc)
        return []


# ---------------------------------------------------------------------------
# GlyTouCan
# ---------------------------------------------------------------------------

async def fetch_glytoucan(glytoucan_id: str) -> dict:
    """Return glycan metadata from GlyTouCan.

    Args:
        glytoucan_id: GlyTouCan accession, e.g. "G00029MO".

    Returns:
        Dict with keys: accession, glycan_name, molecular_formula.
        Empty dict on any error.
    """
    try:
        url = f"https://api.glytoucan.org/core/accession/{glytoucan_id}"
        async with httpx.AsyncClient(timeout=20.0) as client:
            r = await client.get(url)
            r.raise_for_status()
            data = r.json()
        return {
            "accession":        data.get("accession", glytoucan_id),
            "glycan_name":      data.get("name") or data.get("glycan_name", ""),
            "molecular_formula": data.get("molecular_formula") or data.get("formula", ""),
        }
    except Exception as exc:
        _log.debug("fetch_glytoucan failed for %s: %s", glytoucan_id, exc)
        return {}


# ---------------------------------------------------------------------------
# Rule-based N-glycosylation prediction (NetNGlyc sequon rule)
# ---------------------------------------------------------------------------

def _predict_sequon_nglycosylation(sequence: str) -> list[dict]:
    """Rule-based N-glycosylation site prediction using the NxS/T sequon.

    Scans for the canonical N-glycosylation sequon N-x-S/T where x is any
    amino acid except proline (P).

    # NxS/T N-glycosylation sequon: Apweiler R 1999 Biochim Biophys Acta
    # doi:10.1016/S0304-4165(99)00165-8

    Args:
        sequence: Single-letter amino acid sequence string.

    Returns:
        List of dicts: {position (1-based), sequon, predicted_score, method}.
    """
    sites: list[dict] = []
    seq = sequence.upper()
    # Pattern: N followed by any non-P residue followed by S or T
    for match in re.finditer(r"N(?=[^P][ST])", seq):
        pos = match.start()  # 0-based
        sequon = seq[pos: pos + 3] if pos + 3 <= len(seq) else seq[pos:]
        sites.append({
            "position":        pos + 1,          # 1-based
            "sequon":          sequon,
            # Default score for a pure rule-based hit; NetNGlyc threshold is 0.5
            # Citation: Gupta R 2004 Glycobiology doi:10.1093/glycob/cwh148
            "predicted_score": 0.7,
            "method":          "sequon_rule",
        })
    return sites


async def _predict_dtu_nglycosylation(sequence: str) -> list[dict]:
    """Attempt NetNGlyc prediction via the DTU web server.

    Posts to the DTU NetNGlyc 1.0 server and parses the HTML response for
    threshold-0.5 predictions.  Returns an empty list if the server is
    unreachable or returns an unexpected response.

    Args:
        sequence: Single-letter amino acid sequence string.

    Returns:
        List of dicts: {position, sequon, predicted_score, method}.
    """
    try:
        fasta = f">seq\n{sequence}"
        form_data = {
            "seqform":  "FASTA",
            "SEQ":      fasta,
            "orgtype":  "euk",
            "outform":  "-h",
        }
        url = "https://services.healthtech.dtu.dk/cgi-bin/webface2.fcgi"
        async with httpx.AsyncClient(timeout=30.0) as client:
            r = await client.post(url, data=form_data)
            r.raise_for_status()
        html = r.text

        sites: list[dict] = []
        # DTU HTML output rows typically look like:
        #   Asn-X-Ser/Thr sequon  Pos  Potential  Jury  NGlyc+/-
        # We search for lines with a score >= 0.5 (the prediction threshold).
        # Pattern captures: position, sequon-text, score
        row_re = re.compile(
            r"(\d+)\s+[A-Z]\s+[A-Z]\s+[A-Z]\s+(0\.\d+)\s+\+",
            re.IGNORECASE,
        )
        for match in row_re.finditer(html):
            pos   = int(match.group(1))
            score = float(match.group(2))
            # Re-derive sequon from sequence for consistency
            idx  = pos - 1
            sequon = sequence[idx: idx + 3].upper() if idx + 3 <= len(sequence) else sequence[idx:].upper()
            sites.append({
                "position":        pos,
                "sequon":          sequon,
                "predicted_score": score,
                "method":          "NetNGlyc-DTU",
            })
        return sites
    except Exception as exc:
        _log.debug("DTU NetNGlyc server call failed: %s", exc)
        return []


async def predict_nglycosylation(sequence: str) -> list[dict]:
    """Predict N-glycosylation sites, preferring DTU server with sequon fallback.

    Attempts the DTU NetNGlyc server first.  If the server returns no results
    (due to timeout, outage, or parsing failure), falls back to the pure
    rule-based NxS/T sequon scan.

    # NxS/T N-glycosylation sequon: Apweiler R 1999 Biochim Biophys Acta
    # doi:10.1016/S0304-4165(99)00165-8

    Args:
        sequence: Single-letter amino acid sequence string.

    Returns:
        List of prediction dicts; deduplication by position is applied.
    """
    if not sequence:
        return []

    # Always compute rule-based predictions as baseline
    rule_sites   = _predict_sequon_nglycosylation(sequence)
    server_sites = await _predict_dtu_nglycosylation(sequence)

    if server_sites:
        # Merge: use server results; add any rule-based positions not covered
        server_positions = {s["position"] for s in server_sites}
        extras = [s for s in rule_sites if s["position"] not in server_positions]
        return server_sites + extras

    return rule_sites


# ---------------------------------------------------------------------------
# Top-level data aggregator
# ---------------------------------------------------------------------------

async def fetch_glycan_data(uniprot_acc: str, sequence: str) -> dict:
    """Fetch and aggregate glycan data for a protein.

    Calls predict_nglycosylation and fetch_glyconnect concurrently, then
    returns a combined dict for consumption by the glycan analyzer.

    Args:
        uniprot_acc: UniProt accession string.
        sequence:    Single-letter amino acid sequence.

    Returns:
        Dict with keys:
            "predicted_sites"     — list of dicts from predict_nglycosylation
            "known_compositions"  — list of dicts from fetch_glyconnect
    """
    predicted, known = await asyncio.gather(
        predict_nglycosylation(sequence),
        fetch_glyconnect(uniprot_acc),
        return_exceptions=True,
    )

    return {
        "predicted_sites":    predicted if isinstance(predicted, list) else [],
        "known_compositions": known     if isinstance(known, list)     else [],
    }
