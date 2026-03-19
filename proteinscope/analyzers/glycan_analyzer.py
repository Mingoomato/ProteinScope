"""Glycan Analyzer — N- and O-linked glycosylation site analysis.

Integrates experimental glycan data from GlyConnect, rule-based N-glycosylation
predictions (NxS/T sequon), and UniProt glycosylation feature annotations into a
unified GlycanAnalysis report.  Includes host-expression compatibility assessment
and Gemini-powered glycan engineering opportunity synthesis.

Background:
  N-linked glycosylation occurs co-translationally at Asn residues in the
  sequon N-x-S/T (x ≠ P).  O-linked glycosylation is added post-translationally
  to Ser/Thr residues.  Glycan complexity and microheterogeneity are major
  quality-attribute concerns in biopharmaceutical manufacturing.

  # NxS/T sequon: Apweiler R 1999 Biochim Biophys Acta
  # doi:10.1016/S0304-4165(99)00165-8
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

class GlycanSite(BaseModel):
    """A single glycosylation site on the protein."""

    position: int
    sequon: str                                             # e.g. "NCS"
    glycan_type: str                                        # "N-linked" | "O-linked"
    predicted_occupancy: Optional[float] = None            # NetNGlyc score (0–1)
    known_compositions: List[str] = Field(default_factory=list)   # from GlyConnect
    glytoucan_id: Optional[str] = None
    provenance: Optional[DataProvenance] = None


class GlycanAnalysis(BaseModel):
    """Full glycan analysis report for a gene/protein."""

    gene: str
    sites: List[GlycanSite]
    total_n_linked: int
    total_o_linked: int
    host_compatibility_impact: str      # host expression recommendation
    engineering_opportunities: List[str]
    gemini_interpretation: str = ""
    timestamp: datetime


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_uniprot_glycan_features(entry_features: list) -> list[dict]:
    """Extract glycosylation annotations from a UniProt features list.

    Iterates the raw JSON features list returned by the UniProt REST API and
    selects items where featureType == "Glycosylation".  Determines glycan
    linkage type from the description text.

    Args:
        entry_features: Raw list of feature dicts from UniProt entry JSON.

    Returns:
        List of dicts: {position, glycan_type, description, source}.
    """
    sites: list[dict] = []
    for feat in entry_features or []:
        try:
            if str(feat.get("type") or feat.get("featureType", "")).lower() != "glycosylation":
                continue
            desc = str(feat.get("description", "")).lower()
            pos_data = feat.get("location", {}).get("start", {})
            pos = int(pos_data.get("value", 0)) if isinstance(pos_data, dict) else 0
            if pos == 0:
                # Fallback: some responses use "begin" key
                pos = int(feat.get("begin", 0))
            if pos == 0:
                continue
            glycan_type = "N-linked" if "n-linked" in desc else "O-linked"
            sites.append({
                "position":    pos,
                "glycan_type": glycan_type,
                "description": feat.get("description", ""),
                "source":      "UniProt",
            })
        except Exception as exc:
            _log.debug("Could not parse glycan feature: %s — %s", feat, exc)
    return sites


def _assess_host_compatibility(n_linked_count: int, o_linked_count: int) -> str:
    """Derive an expression host recommendation from glycan site counts.

    Rules are based on the requirement for mammalian-type glycosylation
    machinery to produce complex N-glycans and the known O-glycosylation
    patterns of common production hosts.

    # Mammalian expression for complex N-glycans: Jayapal KP 2007
    # Chem Eng Prog doi:10.1002/cep.20070103

    Args:
        n_linked_count: Number of N-linked glycosylation sites.
        o_linked_count: Number of O-linked glycosylation sites.

    Returns:
        A human-readable host compatibility recommendation string.
    """
    parts: list[str] = []

    if n_linked_count > 3:
        # Complex multi-site N-glycosylation requires mammalian machinery
        # Citation: Jayapal KP 2007 Chem Eng Prog doi:10.1002/cep.20070103
        parts.append(
            "CHO-K1 or HEK293 preferred (mammalian glycosylation machinery "
            "required for complex N-glycan processing)"
        )
    elif n_linked_count > 0:
        parts.append(
            "Mammalian (CHO/HEK293) preferred for accurate N-glycan structure; "
            "Pichia pastoris viable with high-mannose glycoforms"
        )
    else:
        # No N-linked sites: prokaryotic or methylotrophic yeast expression feasible
        parts.append("E. coli or Pichia pastoris viable (no N-linked glycosylation required)")

    if o_linked_count > 0:
        # Pichia pastoris adds O-glycans but uses mannose rather than GalNAc
        # Citation: Hamilton SR 2006 Science doi:10.1126/science.1130256
        parts.append(
            "O-linked sites detected: Pichia pastoris can O-glycosylate but "
            "with mannose-based O-glycans — may differ from native GalNAc pattern"
        )

    return "; ".join(parts) if parts else "No expression host constraint identified."


# ---------------------------------------------------------------------------
# Gemini synthesis
# ---------------------------------------------------------------------------

async def _synthesize_glycan(
    gene: str,
    sites: list[GlycanSite],
    host_recommendation: str,
) -> dict:
    """Call Gemini to interpret glycan biology and propose engineering strategies.

    Returns dict with keys: interpretation, engineering_opportunities.
    Returns empty dict on failure.
    """
    try:
        from core.gemini_interpreter import _call

        sites_summary = _json.dumps(
            [
                {
                    "position":           s.position,
                    "sequon":             s.sequon,
                    "glycan_type":        s.glycan_type,
                    "predicted_occupancy": s.predicted_occupancy,
                    "known_compositions": s.known_compositions,
                }
                for s in sites[:30]
            ],
            indent=2,
        )

        prompt = (
            f"You are a glycobiology and protein engineering expert. "
            f"Analyze the following glycosylation sites for {gene}:\n"
            f"{sites_summary}\n\n"
            f"Host expression recommendation: {host_recommendation}\n\n"
            "Provide:\n"
            "1. A concise biological interpretation of the glycan landscape "
            "(significance for folding, stability, receptor binding, immunogenicity).\n"
            "2. A list of specific glycan engineering opportunities — include:\n"
            "   - N-glycan site removal via N→Q mutagenesis (to reduce heterogeneity)\n"
            "   - N-glycan site addition (consensus sequon insertion for stability)\n"
            "   - O-glycan engineering if relevant\n"
            "   - Aglycosylation strategies if appropriate\n\n"
            "Return JSON only:\n"
            "{\n"
            '  "interpretation": "...",\n'
            '  "engineering_opportunities": ["...", "..."]\n'
            "}\n"
            "Return ONLY raw JSON — no markdown fences, no explanation."
        )

        raw = await _call(prompt)
        if not raw:
            return {}

        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        return _json.loads(cleaned)
    except Exception as exc:
        _log.debug("Gemini glycan synthesis failed for %s: %s", gene, exc)
        return {}


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_glycan_analysis(
    gene: str,
    sequence: str,
    entry_features: list,       # UniProt features list (raw JSON)
    glycan_data: dict,          # from fetchers.glycan.fetch_glycan_data()
    step_cb=None,
) -> Optional[GlycanAnalysis]:
    """Run full glycan analysis for a protein.

    Args:
        gene:           Gene symbol, e.g. "ERBB2".
        sequence:       Single-letter amino acid sequence.
        entry_features: Raw UniProt features list for the canonical isoform.
        glycan_data:    Output of fetchers.glycan.fetch_glycan_data().
        step_cb:        Optional async progress callback (receives str).

    Returns:
        GlycanAnalysis report, or None if a fatal error occurs.
    """
    async def _step(msg: str) -> None:
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    try:
        # ── 1. Parse UniProt glycosylation features ──────────────────────────
        await _step("[1/5] Parsing UniProt glycosylation features...")
        uniprot_sites = _parse_uniprot_glycan_features(entry_features)

        # ── 2. Merge UniProt experimental + predicted N-glycosylation sites ──
        await _step("[2/5] Processing predicted N-glycosylation sites...")
        # Build a position-keyed dict; UniProt experimental sites take priority
        merged: dict[int, GlycanSite] = {}

        for us in uniprot_sites:
            pos = us["position"]
            # Derive sequon from sequence (1-based position)
            idx   = pos - 1
            sequon = sequence[idx: idx + 3].upper() if 0 <= idx < len(sequence) else ""
            merged[pos] = GlycanSite(
                position=pos,
                sequon=sequon,
                glycan_type=us["glycan_type"],
                provenance=DataProvenance(
                    source="UniProt",
                    evidence_grade=EvidenceGrade.EXPERIMENTAL,
                    scientific_caveat="Glycosylation site from UniProt annotation.",
                ),
            )

        # Add predicted N-glycosylation sites not already covered by UniProt
        for pred in glycan_data.get("predicted_sites", []):
            pos = int(pred.get("position", 0))
            if pos == 0:
                continue
            if pos not in merged:
                sequon = str(pred.get("sequon", ""))
                if not sequon:
                    idx = pos - 1
                    sequon = sequence[idx: idx + 3].upper() if 0 <= idx < len(sequence) else ""
                method = pred.get("method", "sequon_rule")
                score  = pred.get("predicted_score")
                merged[pos] = GlycanSite(
                    position=pos,
                    sequon=sequon,
                    glycan_type="N-linked",
                    predicted_occupancy=float(score) if score is not None else None,
                    provenance=DataProvenance(
                        source=method,
                        evidence_grade=EvidenceGrade.COMPUTATIONAL,
                        # NetNGlyc threshold 0.5: Gupta R 2004 Glycobiology doi:10.1093/glycob/cwh148
                        scientific_caveat=(
                            "Predicted N-glycosylation site (NetNGlyc/sequon rule); "
                            "experimental validation required."
                        ),
                        method=method,
                    ),
                )

        # ── 3. Match known GlyConnect compositions ───────────────────────────
        await _step("[3/5] Matching known GlyConnect compositions...")
        known_by_pos: dict[int, list[dict]] = {}
        for comp in glycan_data.get("known_compositions", []):
            p = comp.get("position")
            if p is None:
                continue
            try:
                p = int(p)
            except (TypeError, ValueError):
                continue
            known_by_pos.setdefault(p, []).append(comp)

        for pos, site in merged.items():
            if pos in known_by_pos:
                compositions  = [
                    c.get("glycan_composition", "")
                    for c in known_by_pos[pos]
                    if c.get("glycan_composition")
                ]
                glytoucan_ids = [
                    c.get("glytoucan_id", "")
                    for c in known_by_pos[pos]
                    if c.get("glytoucan_id")
                ]
                site.known_compositions = list(dict.fromkeys(compositions))   # dedup, order-preserving
                site.glytoucan_id       = glytoucan_ids[0] if glytoucan_ids else None
                # Upgrade provenance to experimental when GlyConnect data exists
                site.provenance = DataProvenance(
                    source="GlyConnect",
                    evidence_grade=EvidenceGrade.EXPERIMENTAL,
                    scientific_caveat="Glycan composition from GlyConnect experimental database.",
                )

        # Sorted list of sites by position
        sites_sorted = sorted(merged.values(), key=lambda s: s.position)

        # ── 4. Assess host compatibility impact ───────────────────────────────
        await _step("[4/5] Assessing host compatibility impact...")
        total_n = sum(1 for s in sites_sorted if s.glycan_type == "N-linked")
        total_o = sum(1 for s in sites_sorted if s.glycan_type == "O-linked")
        host_rec = _assess_host_compatibility(total_n, total_o)

        # ── 5. Gemini synthesis ───────────────────────────────────────────────
        await _step("[5/5] Gemini synthesis...")
        gemini_data = await _synthesize_glycan(gene, sites_sorted, host_rec)

        interpretation       = str(gemini_data.get("interpretation", "")).strip()
        engineering_opps_raw = gemini_data.get("engineering_opportunities", [])
        engineering_opps: list[str] = (
            [str(o).strip() for o in engineering_opps_raw if o]
            if isinstance(engineering_opps_raw, list)
            else []
        )

        # Provide fallback engineering opportunities if Gemini unavailable
        if not engineering_opps:
            engineering_opps = _default_engineering_opportunities(sites_sorted)

        return GlycanAnalysis(
            gene=gene,
            sites=sites_sorted,
            total_n_linked=total_n,
            total_o_linked=total_o,
            host_compatibility_impact=host_rec,
            engineering_opportunities=engineering_opps,
            gemini_interpretation=interpretation,
            timestamp=datetime.utcnow(),
        )

    except Exception as exc:
        _log.debug("run_glycan_analysis failed for %s: %s", gene, exc)
        return None


def _default_engineering_opportunities(sites: list[GlycanSite]) -> list[str]:
    """Produce rule-based engineering suggestions when Gemini is unavailable."""
    opps: list[str] = []
    n_sites = [s for s in sites if s.glycan_type == "N-linked"]
    o_sites = [s for s in sites if s.glycan_type == "O-linked"]

    if len(n_sites) > 1:
        positions_str = ", ".join(f"N{s.position}" for s in n_sites[:3])
        opps.append(
            f"N→Q mutagenesis at {positions_str} to reduce N-glycan "
            "microheterogeneity for analytical homogeneity."
        )
    if n_sites:
        opps.append(
            "Consider site-directed introduction of additional N-glycosylation "
            "sequons (N-x-S/T) to improve thermal stability or mask immunogenic epitopes."
        )
    if o_sites:
        opps.append(
            "O-linked sites identified; consider Ser/Thr→Ala mutagenesis to "
            "ablate O-glycosylation if mammalian host O-glycan patterns are problematic."
        )
    if not opps:
        opps.append(
            "No glycan engineering interventions immediately indicated; "
            "confirm aglycosylation compatibility for chosen expression host."
        )
    return opps
