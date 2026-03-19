"""SynBioHub iGEM Registry parts fetcher.

Searches for BioBrick parts by function keyword via the SynBioHub public REST API.
REST: GET https://synbiohub.org/public/igem/search?q={keyword}&objectType=ComponentDefinition
Returns: [{part_id, name, description, uri, part_type}]
Fallback: hardcoded canonical iGEM BioBrick list by keyword category
"""

from __future__ import annotations

import httpx


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

async def search_igem_parts(function_keyword: str) -> list[dict]:
    """Search SynBioHub/iGEM Registry for BioBrick ComponentDefinitions.

    Args:
        function_keyword: Free-text keyword to search (e.g. "promoter", "RBS",
                          "terminator", "repressor", "fluorescent reporter").

    Returns:
        List of dicts with keys: part_id, name, description, uri, part_type.
        Falls back to ``_get_default_parts`` on any network or parse error.
    """
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            r = await client.get(
                "https://synbiohub.org/public/igem/search",
                params={
                    "q": function_keyword,
                    "objectType": "ComponentDefinition",
                    "limit": 10,
                },
                headers={"Accept": "application/json"},
            )
            r.raise_for_status()
            raw = r.json()

            parts: list[dict] = []
            # SynBioHub returns a list of objects directly or wrapped under a key
            items = raw if isinstance(raw, list) else raw.get("results", [])
            for item in items:
                if not isinstance(item, dict):
                    continue
                uri = str(item.get("uri", "") or "")
                display_id = str(item.get("displayId", "") or "")
                name = str(item.get("name", "") or display_id or "")
                description = str(item.get("description", "") or "")
                # Infer part_id: prefer displayId; fall back to last URI segment
                if display_id:
                    part_id = display_id
                elif uri:
                    part_id = uri.rstrip("/").split("/")[-1]
                else:
                    continue
                # Infer type from the URI path or displayId prefix
                part_type = _infer_part_type(uri, display_id, function_keyword)
                parts.append({
                    "part_id": part_id,
                    "name": name,
                    "description": description,
                    "uri": uri,
                    "part_type": part_type,
                })

            # If the live query returned nothing useful, fall back to defaults
            if not parts:
                return _get_default_parts(function_keyword)
            return parts

    except Exception:
        return _get_default_parts(function_keyword)


# ---------------------------------------------------------------------------
# Type inference helper
# ---------------------------------------------------------------------------

def _infer_part_type(uri: str, display_id: str, keyword: str) -> str:
    """Heuristically infer the BioBrick part type from URI / displayId."""
    combined = (uri + display_id + keyword).lower()
    if "promoter" in combined:
        return "Promoter"
    if "terminator" in combined:
        return "Terminator"
    if "rbs" in combined or "ribosome" in combined:
        return "RBS"
    if "reporter" in combined or "gfp" in combined or "rfp" in combined or "mcherry" in combined:
        return "Reporter"
    if "cds" in combined or "coding" in combined:
        return "CDS"
    if "regulatory" in combined or "repressor" in combined or "activator" in combined:
        return "Regulatory"
    return "ComponentDefinition"


# ---------------------------------------------------------------------------
# Hardcoded fallback library
# ---------------------------------------------------------------------------

def _get_default_parts(keyword: str) -> list[dict]:
    """Return a curated list of well-characterised iGEM BioBrick parts.

    Used when the SynBioHub REST call fails or returns no results.
    Covers the most common functional categories in synthetic biology.

    References:
        iGEM Parts Registry: http://parts.igem.org
        SynBioHub: https://synbiohub.org/public/igem
    """
    kw = keyword.lower()

    # ── Promoters ────────────────────────────────────────────────────────────
    if "promoter" in kw:
        return [
            {
                "part_id": "BBa_R0040",
                "name": "Ptrc",
                "description": "IPTG-inducible trc promoter; hybrid trp/lac; strong inducible expression in E. coli",
                "uri": "https://synbiohub.org/public/igem/BBa_R0040",
                "part_type": "Promoter",
            },
            {
                "part_id": "BBa_J23100",
                "name": "Constitutive promoter (Anderson series, strongest)",
                "description": "Anderson constitutive promoter family member; highest strength in the J231xx series",
                "uri": "https://synbiohub.org/public/igem/BBa_J23100",
                "part_type": "Promoter",
            },
            {
                "part_id": "BBa_J23106",
                "name": "Constitutive promoter (Anderson series, medium-high)",
                "description": "Anderson constitutive promoter; ~68% relative strength vs J23100; E. coli",
                "uri": "https://synbiohub.org/public/igem/BBa_J23106",
                "part_type": "Promoter",
            },
            {
                "part_id": "BBa_I14018",
                "name": "Plac",
                "description": "LacI-repressed promoter; derepressed by IPTG; widely used lac-derived promoter",
                "uri": "https://synbiohub.org/public/igem/BBa_I14018",
                "part_type": "Promoter",
            },
            {
                "part_id": "BBa_R0011",
                "name": "Plac/ara-1",
                "description": "Hybrid lac/ara promoter; induced by arabinose; tighter repression than Plac",
                "uri": "https://synbiohub.org/public/igem/BBa_R0011",
                "part_type": "Promoter",
            },
        ]

    # ── Terminators ──────────────────────────────────────────────────────────
    if "terminator" in kw:
        return [
            {
                "part_id": "BBa_B0015",
                "name": "Double terminator",
                "description": "Strong bidirectional double terminator (rrnB T1-T2); gold standard for E. coli constructs",
                "uri": "https://synbiohub.org/public/igem/BBa_B0015",
                "part_type": "Terminator",
            },
            {
                "part_id": "BBa_B0010",
                "name": "rrnB T1 terminator",
                "description": "rrnB operon T1 terminator; high efficiency intrinsic terminator in E. coli",
                "uri": "https://synbiohub.org/public/igem/BBa_B0010",
                "part_type": "Terminator",
            },
            {
                "part_id": "BBa_B0012",
                "name": "rrnB T2 terminator",
                "description": "rrnB operon T2 terminator; often used downstream of T1 for complete transcriptional stop",
                "uri": "https://synbiohub.org/public/igem/BBa_B0012",
                "part_type": "Terminator",
            },
            {
                "part_id": "BBa_B0054",
                "name": "Transcriptional terminator",
                "description": "Strong synthetic terminator; characterised in iGEM 2006",
                "uri": "https://synbiohub.org/public/igem/BBa_B0054",
                "part_type": "Terminator",
            },
        ]

    # ── Ribosome Binding Sites ────────────────────────────────────────────────
    if "rbs" in kw or "ribosome" in kw:
        return [
            {
                "part_id": "BBa_B0034",
                "name": "Medium strength RBS (Elowitz)",
                "description": "Most commonly used RBS in iGEM constructs; ~1000 AU expression in standard conditions",
                "uri": "https://synbiohub.org/public/igem/BBa_B0034",
                "part_type": "RBS",
            },
            {
                "part_id": "BBa_B0030",
                "name": "Strong RBS",
                "description": "Strong RBS derived from T7 gene 10; high translation initiation rate",
                "uri": "https://synbiohub.org/public/igem/BBa_B0030",
                "part_type": "RBS",
            },
            {
                "part_id": "BBa_B0031",
                "name": "Weak RBS",
                "description": "Weak Shine-Dalgarno sequence; useful for tuning expression to low levels",
                "uri": "https://synbiohub.org/public/igem/BBa_B0031",
                "part_type": "RBS",
            },
            {
                "part_id": "BBa_B0032",
                "name": "Medium RBS",
                "description": "Medium-strength RBS; ~30% relative strength vs BBa_B0034",
                "uri": "https://synbiohub.org/public/igem/BBa_B0032",
                "part_type": "RBS",
            },
        ]

    # ── Reporters ────────────────────────────────────────────────────────────
    if "reporter" in kw or "gfp" in kw or "fluorescent" in kw:
        return [
            {
                "part_id": "BBa_E0040",
                "name": "GFP (green fluorescent protein)",
                "description": "Enhanced GFP coding sequence; excitation 395/475 nm, emission 509 nm",
                "uri": "https://synbiohub.org/public/igem/BBa_E0040",
                "part_type": "Reporter",
            },
            {
                "part_id": "BBa_E1010",
                "name": "mRFP1 (monomeric red fluorescent protein)",
                "description": "Monomeric RFP1; excitation 584 nm, emission 607 nm; widely used for dual-reporter assays",
                "uri": "https://synbiohub.org/public/igem/BBa_E1010",
                "part_type": "Reporter",
            },
            {
                "part_id": "BBa_J04450",
                "name": "mCherry",
                "description": "Bright monomeric red fluorescent protein; photostable; excitation 587 nm, emission 610 nm",
                "uri": "https://synbiohub.org/public/igem/BBa_J04450",
                "part_type": "Reporter",
            },
        ]

    # ── Regulatory / toggle / logic ───────────────────────────────────────────
    if "regulatory" in kw or "repressor" in kw or "toggle" in kw or "gate" in kw:
        return [
            {
                "part_id": "BBa_C0012",
                "name": "LacI repressor",
                "description": "LacI coding sequence; represses Plac/Ptrc; derepressed by IPTG",
                "uri": "https://synbiohub.org/public/igem/BBa_C0012",
                "part_type": "Regulatory",
            },
            {
                "part_id": "BBa_C0051",
                "name": "cI repressor (lambda phage)",
                "description": "Lambda phage cI repressor; represses PR/PL; used in toggle switches",
                "uri": "https://synbiohub.org/public/igem/BBa_C0051",
                "part_type": "Regulatory",
            },
            {
                "part_id": "BBa_C0080",
                "name": "TetR repressor",
                "description": "Tetracycline repressor; represses Ptet; derepressed by aTc (anhydrotetracycline)",
                "uri": "https://synbiohub.org/public/igem/BBa_C0080",
                "part_type": "Regulatory",
            },
        ]

    # ── Default: return empty list (caller must handle) ──────────────────────
    return []
