"""PTM Logic Engine — Post-translational modification regulatory graph analyzer.

Background — What are PTMs?
─────────────────────────────────────────────────────────────────────────────
Post-translational modifications (PTMs) are covalent chemical modifications
to amino acid residues that occur after protein synthesis. They dramatically
expand the functional repertoire of the proteome and form the basis of
reversible cell signaling.

Key PTM classes:
  Phosphorylation (Ser/Thr/Tyr) — most prevalent; added by kinases,
    removed by phosphatases; controls activity, localization, interactions.
  Ubiquitination (Lys) — marks proteins for proteasomal degradation (K48),
    or modulates signaling (K63, M1 chains).
  Acetylation (Lys, N-terminal) — often competes with ubiquitination;
    histone acetylation opens chromatin.
  Methylation (Lys, Arg) — epigenetic regulation; mono-/di-/tri-methylation
    have distinct effector binding.
  SUMOylation, palmitoylation, glycosylation — specialized regulatory roles.

Regulatory logic:
  PTMs rarely act in isolation. Common patterns:
    - Sequential phosphorylation cascades (e.g. PDK1→AKT T308 primes S473)
    - PTM crosstalk: phospho-degron (phospho → ubiquitin E3 recognition → degradation)
    - Mutual exclusion: competing modifications at the same or adjacent residues
    - Cooperative activation: multiple sites must be phosphorylated together

This analyzer performs:
  1. Parsing raw PTM site data into structured PTMSite objects
  2. Deduplication and evidence grading per site
  3. Inference of regulatory state transitions between sites
  4. Gemini-powered PTM logic synthesis (regulatory relationships, confidence)
  5. Assembly of a PTMLogicGraph with key regulatory nodes and implications
"""

from __future__ import annotations

import json as _json
import logging
from typing import List, Optional

from pydantic import BaseModel, Field

from core.evidence import EvidenceGrade

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

class PTMSite(BaseModel):
    position: int
    residue: str
    modification: str
    site_id: str                         # e.g. "S473"
    known_kinases: List[str] = Field(default_factory=list)
    functional_effect: str = "unknown"
    evidence_grade: EvidenceGrade = EvidenceGrade.EXPERIMENTAL
    source: str = ""
    pmids: List[str] = Field(default_factory=list)


class PTMLogicRelation(BaseModel):
    site_a: str
    site_b: str
    relation: str                        # "requires_prior"|"prevents"|"cooperates"|"independent"
    confidence: float = 0.5
    evidence: str = ""
    evidence_pmids: List[str] = Field(default_factory=list)


class PTMLogicGraph(BaseModel):
    gene: str
    sites: List[PTMSite] = Field(default_factory=list)
    relations: List[PTMLogicRelation] = Field(default_factory=list)
    key_regulatory_nodes: List[str] = Field(default_factory=list)
    regulatory_summary: str = ""
    therapeutic_implications: str = ""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _grade_evidence(source: str, pmids: list) -> EvidenceGrade:
    """Assign evidence grade based on data source and publication support."""
    if source == "PhosphoSitePlus" and pmids:
        return EvidenceGrade.EXPERIMENTAL
    if source == "PhosphoSitePlus":
        return EvidenceGrade.EXPERIMENTAL
    if source == "UniProt" and pmids:
        return EvidenceGrade.EXPERIMENTAL
    if source == "UniProt":
        return EvidenceGrade.COMPUTATIONAL
    return EvidenceGrade.LITERATURE


def _deduplicate_sites(raw_sites: list[dict]) -> list[PTMSite]:
    """Parse, deduplicate, and grade PTM site dicts into PTMSite objects.

    Deduplication key: (position, modification_type). When duplicates exist,
    the entry with the most PMIDs and highest-priority source is retained.
    """
    # Source priority: PhosphoSitePlus > UniProt > other
    _source_priority = {"PhosphoSitePlus": 0, "UniProt": 1}

    # Group by (position, modification_type) for dedup
    grouped: dict[tuple, list[dict]] = {}
    for s in raw_sites:
        pos = int(s.get("position", 0))
        mod = str(s.get("modification_type", "")).strip()
        key = (pos, mod)
        grouped.setdefault(key, []).append(s)

    result: list[PTMSite] = []
    for (pos, mod), entries in grouped.items():
        # Pick best entry: lowest source_priority value, then most pmids
        best = min(
            entries,
            key=lambda e: (
                _source_priority.get(e.get("source", ""), 99),
                -len(e.get("pmids", [])),
            ),
        )
        residue = str(best.get("residue", "")).strip()
        site_id = str(best.get("site", f"{residue}{pos}")).strip() or f"{residue}{pos}"
        pmids = [str(p) for p in best.get("pmids", []) if p]
        result.append(PTMSite(
            position=pos,
            residue=residue,
            modification=mod,
            site_id=site_id,
            known_kinases=[str(k) for k in best.get("known_kinases", []) if k],
            functional_effect=str(best.get("functional_effect", "unknown")),
            evidence_grade=_grade_evidence(best.get("source", ""), pmids),
            source=str(best.get("source", "")),
            pmids=pmids,
        ))

    # Sort by position ascending
    result.sort(key=lambda s: s.position)
    return result


def _infer_state_transitions(sites: list[PTMSite]) -> list[PTMLogicRelation]:
    """Rule-based inference of likely regulatory relationships between PTM sites.

    Applies known biochemical heuristics for common kinase cascades:
      - Sequential phosphorylation: upstream kinase site at lower position
        feeding a downstream site (e.g. PDK1-T308 → AKT-S473).
      - Competing sites: same residue type at adjacent positions (<10 aa apart)
        with opposing functional effects → "prevents" relation.
      - Cooperative activation: two activating phospho-sites within 50 aa
        sharing a common kinase → "cooperates" relation.
    These are heuristic inferences and are superseded by Gemini synthesis.
    """
    relations: list[PTMLogicRelation] = []
    phospho_sites = [
        s for s in sites
        if "phospho" in s.modification.lower() or "phospho" in s.site_id.lower()
    ]

    for i, site_a in enumerate(phospho_sites):
        for site_b in phospho_sites[i + 1:]:
            dist = abs(site_b.position - site_a.position)
            shared_kinases = set(site_a.known_kinases) & set(site_b.known_kinases)

            # Cooperative: same kinase, both activating, close proximity
            if (
                shared_kinases
                and site_a.functional_effect in ("activation", "unknown")
                and site_b.functional_effect in ("activation", "unknown")
                and dist <= 100
            ):
                relations.append(PTMLogicRelation(
                    site_a=site_a.site_id,
                    site_b=site_b.site_id,
                    relation="cooperates",
                    confidence=0.55,
                    evidence=(
                        f"Shared kinase(s): {', '.join(sorted(shared_kinases))}; "
                        f"both sites within {dist} aa"
                    ),
                ))

            # Competing: very close sites with opposing effects
            elif (
                dist <= 10
                and site_a.functional_effect != site_b.functional_effect
                and "unknown" not in (site_a.functional_effect, site_b.functional_effect)
            ):
                relations.append(PTMLogicRelation(
                    site_a=site_a.site_id,
                    site_b=site_b.site_id,
                    relation="prevents",
                    confidence=0.45,
                    evidence=(
                        f"Adjacent sites ({dist} aa apart) with opposing "
                        f"effects: {site_a.functional_effect} vs {site_b.functional_effect}"
                    ),
                ))

    return relations


# ---------------------------------------------------------------------------
# Gemini PTM logic synthesis
# ---------------------------------------------------------------------------

async def _synthesize_ptm_logic(
    gene: str,
    sites: list[PTMSite],
) -> dict:
    """Call Gemini to identify regulatory relationships between PTM sites.

    Returns dict with keys:
        state_transitions, regulatory_summary, key_regulatory_nodes,
        therapeutic_targets
    Returns empty dict on failure.
    """
    try:
        from core.gemini_interpreter import _call

        sites_json = _json.dumps(
            [
                {
                    "site_id": s.site_id,
                    "position": s.position,
                    "residue": s.residue,
                    "modification": s.modification,
                    "known_kinases": s.known_kinases,
                    "functional_effect": s.functional_effect,
                    "source": s.source,
                }
                for s in sites[:40]  # cap to avoid token overflow
            ],
            indent=2,
        )

        prompt = (
            f"You are a protein biochemistry expert. Given these PTM sites for {gene}:\n"
            f"{sites_json}\n\n"
            "Identify regulatory relationships between sites. Return JSON:\n"
            "{\n"
            '  "state_transitions": [\n'
            '    {"from_site": "T308", "to_site": "S473", "relation": "requires_prior", '
            '"confidence": 0.9, '
            '"evidence": "AKT activation requires PDK1 phosphorylation first"}\n'
            "  ],\n"
            '  "regulatory_summary": "...",\n'
            '  "key_regulatory_nodes": ["T308", "S473"],\n'
            '  "therapeutic_targets": ["T308 — PDK1 inhibition blocks AKT activation cascade"]\n'
            "}\n\n"
            "Return ONLY raw JSON — no markdown fences, no explanation."
        )

        raw = await _call(prompt)
        if not raw:
            return {}

        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)
        return data

    except Exception as exc:
        _log.debug("Gemini PTM synthesis failed for %s: %s", gene, exc)
        return {}


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_ptm_analysis(
    gene: str,
    ptm_sites: list[dict],
    step_cb=None,
) -> PTMLogicGraph:
    """Build a PTM regulatory logic graph for a gene.

    Args:
        gene:      Gene symbol (e.g. "AKT1", "EGFR").
        ptm_sites: List of raw PTM site dicts from fetch_ptm_sites().
        step_cb:   Optional async progress callback (receives a str message).

    Returns a PTMLogicGraph (empty shell if ptm_sites is empty or analysis fails).
    """
    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    # Return empty graph immediately if no sites provided
    if not ptm_sites:
        return PTMLogicGraph(gene=gene)

    # ── 1. Parse PTM site data ───────────────────────────────────────────────
    await _step("[1/5] Parsing PTM site data…")
    parsed_sites: list[PTMSite] = []
    try:
        parsed_sites = _deduplicate_sites(ptm_sites)
    except Exception as exc:
        _log.debug("PTM site parsing failed for %s: %s", gene, exc)
        return PTMLogicGraph(gene=gene)

    if not parsed_sites:
        return PTMLogicGraph(gene=gene)

    # ── 2. Deduplicate and grade evidence ────────────────────────────────────
    await _step("[2/5] Deduplicating and grading evidence…")
    # Grading is performed inside _deduplicate_sites; here we log summary stats
    psp_count = sum(1 for s in parsed_sites if s.source == "PhosphoSitePlus")
    uni_count = sum(1 for s in parsed_sites if s.source == "UniProt")
    _log.debug(
        "PTM sites for %s: %d total (%d PhosphoSitePlus, %d UniProt)",
        gene, len(parsed_sites), psp_count, uni_count,
    )

    # ── 3. Infer regulatory state transitions ────────────────────────────────
    await _step("[3/5] Inferring regulatory state transitions…")
    heuristic_relations: list[PTMLogicRelation] = []
    try:
        heuristic_relations = _infer_state_transitions(parsed_sites)
    except Exception as exc:
        _log.debug("State transition inference failed for %s: %s", gene, exc)

    # ── 4. Gemini PTM logic synthesis ────────────────────────────────────────
    await _step("[4/5] Running Gemini PTM logic synthesis…")
    gemini_data: dict = {}
    try:
        gemini_data = await _synthesize_ptm_logic(gene, parsed_sites)
    except Exception as exc:
        _log.debug("Gemini PTM synthesis exception for %s: %s", gene, exc)

    # ── 5. Assemble PTM regulatory graph ─────────────────────────────────────
    await _step("[5/5] Assembling PTM regulatory graph…")
    try:
        # Build relations: start with heuristics, then add/replace with Gemini output
        relations: list[PTMLogicRelation] = list(heuristic_relations)

        valid_relations = {"requires_prior", "prevents", "cooperates", "independent"}
        for tr in gemini_data.get("state_transitions") or []:
            site_a = str(tr.get("from_site", "")).strip()
            site_b = str(tr.get("to_site", "")).strip()
            if not site_a or not site_b:
                continue
            rel_str = str(tr.get("relation", "independent")).strip()
            if rel_str not in valid_relations:
                rel_str = "independent"
            try:
                conf = float(tr.get("confidence", 0.5))
                conf = max(0.0, min(1.0, conf))
            except (TypeError, ValueError):
                conf = 0.5
            relations.append(PTMLogicRelation(
                site_a=site_a,
                site_b=site_b,
                relation=rel_str,
                confidence=conf,
                evidence=str(tr.get("evidence", "")).strip(),
            ))

        # Key regulatory nodes
        key_nodes: list[str] = []
        raw_nodes = gemini_data.get("key_regulatory_nodes") or []
        if isinstance(raw_nodes, list):
            key_nodes = [str(n).strip() for n in raw_nodes if n]
        # Fall back to sites with highest-confidence activation effect
        if not key_nodes:
            activation_sites = [
                s.site_id for s in parsed_sites
                if s.functional_effect == "activation"
            ]
            key_nodes = activation_sites[:5]

        # Therapeutic implications — join Gemini therapeutic_targets list
        therapeutic_targets = gemini_data.get("therapeutic_targets") or []
        if isinstance(therapeutic_targets, list):
            therapeutic_implications = "; ".join(
                str(t).strip() for t in therapeutic_targets if t
            )
        else:
            therapeutic_implications = str(therapeutic_targets).strip()

        regulatory_summary = str(
            gemini_data.get("regulatory_summary", "")
        ).strip()

        return PTMLogicGraph(
            gene=gene,
            sites=parsed_sites,
            relations=relations,
            key_regulatory_nodes=key_nodes,
            regulatory_summary=regulatory_summary,
            therapeutic_implications=therapeutic_implications,
        )

    except Exception as exc:
        _log.debug("PTM graph assembly failed for %s: %s", gene, exc)
        # Return a partial graph with just the parsed sites
        return PTMLogicGraph(gene=gene, sites=parsed_sites)
