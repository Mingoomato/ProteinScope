"""Target Discovery Analyzer — Research Orientation Query (ROQ) handler.

Triggered when a user query describes a disease/pathogen context but no specific
protein. Maps disease → organism (via Gemini), fetches the pathogen proteome,
scores each protein as a drug target, and returns a ranked candidate report.

Pipeline:
    Step 1: Gemini parse → {organism_scientific, disease_name, modality_hint, ...}
    Step 2: fetch_organism_proteome() → raw protein list
    Step 3: For each candidate — compute essentiality score + human homolog check
    Step 4: Rank candidates by target_score
    Step 5: Gemini synthesis → ranked rationale + overall strategy

Returns TargetDiscoveryReport. Never raises — returns minimal report on failure.
"""
from __future__ import annotations

import asyncio
import logging
from datetime import datetime
from typing import Optional, List

import httpx
from pydantic import BaseModel, Field

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Known high-value targets per pathogen genus
# Citation: Odds (2010) Antifungal agents: targets or systems. Mol Microbiol
# Citation: Rex (2016) Antifungal drug resistance. Clin Infect Dis
# ---------------------------------------------------------------------------
_KNOWN_ANTIFUNGAL_TARGETS: dict[str, list[str]] = {
    "Malassezia": ["LIP1", "LIP2", "LIP3", "MGL1", "FAS1", "ERG11", "ERG6", "CYP51"],
    "Candida":    ["ERG11", "FKS1", "CDR1", "MDR1", "SAP1", "SAP2", "HWP1", "ALS3"],
    "Aspergillus":["CYP51A", "FKS1", "GEL1", "CHS1", "PksP", "Lac1"],
    "Trichophyton":["ERG11", "SQS1", "HMG1", "CHS1", "CHS2"],
    "Cryptococcus":["ERG11", "FKS1", "LAC1", "PKA1", "CNA1"],
}

_KNOWN_ANTIBACTERIAL_TARGETS: dict[str, list[str]] = {
    "Staphylococcus": ["MecA", "PBP2a", "GyrA", "GyrB", "FtsZ", "FabI"],
    "Escherichia":    ["GyrA", "GyrB", "FabI", "MurA", "FtsZ", "AcrB"],
    "Pseudomonas":    ["GyrA", "MexAB", "OprM", "AlgD", "FliC", "PilA"],
    "Mycobacterium":  ["InhA", "KatG", "RpoB", "EmbB", "GyrA", "DprE1"],
    "Helicobacter":   ["CagA", "VacA", "HpUreB", "GyrA", "PBP2"],
}

# GO term keywords indicating essentiality / druggability
_ESSENTIAL_GO_KEYWORDS = {
    "dna replication", "rna polymerase", "ribosom", "cell division", "cell wall",
    "membrane synthesis", "fatty acid", "ergosterol", "lipase", "protease",
    "kinase", "atp", "translation", "transcription", "dna repair",
}
_MEMBRANE_GO_KEYWORDS = {
    "membrane", "transmembrane", "cell wall", "cell surface", "secreted",
    "extracellular", "plasma membrane",
}


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class TargetCandidate(BaseModel):
    rank: int
    accession: str
    gene_name: str
    protein_name: str
    organism: str
    target_score: float                           # 0.0–1.0
    essentiality_note: str
    human_identity_pct: Optional[float] = None
    selectivity_risk: bool = False                # True if human_identity_pct > 60%
    known_target: bool = False                    # True if in _KNOWN_*_TARGETS dict
    known_inhibitor_class: Optional[str] = None
    modality_hint: str = "small molecule"
    gemini_rationale: str = ""


class TargetDiscoveryReport(BaseModel):
    query: str
    disease_name: str
    organism_scientific: str
    organism_common: str
    candidates: List[TargetCandidate] = Field(default_factory=list)
    context_summary: str = ""
    gemini_strategy: str = ""
    timestamp: datetime = Field(default_factory=datetime.utcnow)


# ---------------------------------------------------------------------------
# Known inhibitor class lookup
# ---------------------------------------------------------------------------
_INHIBITOR_CLASSES: dict[str, str] = {
    "ERG11": "Azoles (fluconazole, voriconazole)",
    "CYP51": "Azoles (fluconazole, voriconazole)",
    "CYP51A": "Azoles (voriconazole, isavuconazole)",
    "FKS1": "Echinocandins (caspofungin, micafungin)",
    "LIP1": "Lipase inhibitors (orlistat analogs)",
    "LIP2": "Lipase inhibitors (orlistat analogs)",
    "LIP3": "Lipase inhibitors (orlistat analogs)",
    "FAS1": "Fatty acid synthase inhibitors (cerulenin, C75)",
    "ERG6": "Sterol methyltransferase inhibitors",
    "GyrA": "Fluoroquinolones (ciprofloxacin, levofloxacin)",
    "GyrB": "Aminocoumarins (novobiocin)",
    "FabI": "Triclosan, diazaborines",
    "MurA": "Fosfomycin",
    "FtsZ": "Berberine, PC190723 (investigational)",
    "InhA": "Isoniazid (via KatG activation)",
    "RpoB": "Rifamycins (rifampicin)",
    "EmbB": "Ethambutol",
}


# ---------------------------------------------------------------------------
# Human homolog check
# ---------------------------------------------------------------------------

async def _check_human_homolog(gene_name: str, pathogen_sequence: str) -> Optional[float]:
    """Fetch human ortholog sequence from UniProt and compute identity %.

    Returns identity percentage (0-100) or None if not found.
    """
    if not gene_name or not pathogen_sequence:
        return None
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            resp = await client.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={
                    "query": f"gene_exact:{gene_name} AND organism_id:9606 AND reviewed:true",
                    "fields": "accession,sequence",
                    "format": "json",
                    "size": 1,
                },
            )
            if resp.status_code != 200:
                return None
            data = resp.json()
            results = data.get("results") or []
            if not results:
                return None
            human_seq = results[0].get("sequence", {}).get("value", "")
            if not human_seq:
                return None
        # Pairwise alignment
        from analyzers.sequence_aligner import align_pairwise
        result = align_pairwise(pathogen_sequence, human_seq)
        if result is None:
            return None
        return result.identity_pct
    except Exception as exc:
        _log.debug("Human homolog check failed for %s: %s", gene_name, exc)
        return None


async def _fetch_sequence(accession: str) -> str:
    """Fetch FASTA sequence for a UniProt accession."""
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            resp = await client.get(
                f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
            )
            if resp.status_code != 200:
                return ""
            lines = resp.text.splitlines()
            return "".join(l for l in lines if not l.startswith(">"))
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def _score_essentiality(protein: dict) -> tuple[float, str]:
    """Return (essentiality_weight 0.4-1.0, note string).

    # Citation: Odds FC (2010) Antimicrob Agents Chemother 54:2777. doi:10.1128/AAC.01345-09
    """
    go_terms_lower = " ".join(protein.get("go_terms") or []).lower()
    func_lower = (protein.get("function_summary") or "").lower()
    pname_lower = (protein.get("protein_name") or "").lower()
    combined = go_terms_lower + " " + func_lower + " " + pname_lower

    if any(kw in combined for kw in _ESSENTIAL_GO_KEYWORDS):
        if any(kw in combined for kw in _MEMBRANE_GO_KEYWORDS):
            return 1.0, "Essential membrane/cell-wall function"
        return 0.9, "Essential metabolic/replication function"
    if any(kw in combined for kw in _MEMBRANE_GO_KEYWORDS):
        return 0.7, "Membrane-associated (accessible target)"
    if "secreted" in combined or "extracellular" in combined:
        return 0.65, "Secreted / extracellular protein"
    return 0.4, "Function not clearly essential"


def _is_known_target(gene_name: str, organism_scientific: str) -> tuple[bool, Optional[str]]:
    """Check if gene is in known target lists. Returns (is_known, inhibitor_class)."""
    genus = (organism_scientific.split() or [""])[0]
    known_list: list[str] = []
    for genus_key, genes in {**_KNOWN_ANTIFUNGAL_TARGETS, **_KNOWN_ANTIBACTERIAL_TARGETS}.items():
        if genus_key.lower() in genus.lower():
            known_list.extend(genes)

    gene_upper = gene_name.upper()
    if gene_upper in [g.upper() for g in known_list]:
        inhibitor = _INHIBITOR_CLASSES.get(gene_upper)
        return True, inhibitor
    return False, None


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

async def run_target_discovery(
    query: str,
    step_cb=None,
) -> TargetDiscoveryReport:
    """Full target discovery pipeline for Research Orientation Queries.

    Args:
        query:    Raw user query text (language-agnostic)
        step_cb:  Optional async callable(message: str) for SSE progress

    Returns:
        TargetDiscoveryReport — never raises
    """

    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    try:
        # Step 1: Gemini parse
        await _step("Identifying target organism and disease context...")
        from core.gemini_interpreter import parse_target_discovery_query
        parsed = await parse_target_discovery_query(query)

        organism_scientific = parsed.get("organism_scientific", "unknown")
        organism_common = parsed.get("organism_common", "unknown")
        disease_name = parsed.get("disease_name", "unknown")
        modality_hint = parsed.get("modality_hint", "small molecule")
        context_summary = parsed.get("context_summary", "")

        if organism_scientific == "unknown":
            return TargetDiscoveryReport(
                query=query,
                disease_name=disease_name,
                organism_scientific=organism_scientific,
                organism_common=organism_common,
                context_summary="Could not identify target organism from query.",
                gemini_strategy="Please specify the pathogen or disease organism more precisely.",
            )

        # Step 2: Fetch proteome
        await _step(f"Fetching {organism_scientific} proteome from UniProt...")
        from fetchers.uniprot_proteome import fetch_organism_proteome
        raw_proteins = await fetch_organism_proteome(organism_scientific, max_proteins=80)

        if not raw_proteins:
            return TargetDiscoveryReport(
                query=query,
                disease_name=disease_name,
                organism_scientific=organism_scientific,
                organism_common=organism_common,
                context_summary=context_summary,
                gemini_strategy=(
                    f"No UniProt entries found for {organism_scientific}. "
                    "Try using the full scientific name (genus + species)."
                ),
            )

        # Step 3: Score candidates (with concurrency cap for human homolog checks)
        await _step(f"Scoring {len(raw_proteins)} candidates for druggability and selectivity...")
        semaphore = asyncio.Semaphore(3)

        scored: list[dict] = []

        async def _process_one(protein: dict) -> Optional[dict]:
            gene_name = protein.get("gene_name", "")
            accession = protein.get("accession", "")

            # Essentiality
            ess_weight, ess_note = _score_essentiality(protein)

            # Known target bonus
            is_known, inhibitor_class = _is_known_target(gene_name, organism_scientific)
            known_bonus = 0.3 if is_known else 0.0

            # Human homolog check (rate-limited)
            human_identity: Optional[float] = None
            selectivity_risk = False
            async with semaphore:
                seq = await _fetch_sequence(accession)
                if seq and gene_name:
                    human_identity = await _check_human_homolog(gene_name, seq)

            if human_identity is not None:
                if human_identity > 60.0:
                    selectivity_bonus = 0.3
                    selectivity_risk = True
                elif human_identity > 30.0:
                    selectivity_bonus = 1.0
                else:
                    selectivity_bonus = 1.5   # highly selective
            else:
                selectivity_bonus = 1.0  # unknown → neutral

            # Literature score heuristic: function_summary length as proxy
            func_len = len(protein.get("function_summary") or "")
            literature_score = min(1.0, func_len / 300.0)

            # Final score
            # Citation: Russ & Lampel (2005) Drug Discov Today 10:1607. doi:10.1016/S1359-6446(05)03666-4
            raw_score = (ess_weight * selectivity_bonus * (0.5 + 0.5 * literature_score)) + known_bonus
            target_score = round(min(1.0, raw_score), 3)

            return {
                "accession": accession,
                "gene_name": gene_name,
                "protein_name": protein.get("protein_name", ""),
                "organism": organism_scientific,
                "target_score": target_score,
                "essentiality_note": ess_note,
                "human_identity_pct": round(human_identity, 1) if human_identity is not None else None,
                "selectivity_risk": selectivity_risk,
                "known_target": is_known,
                "known_inhibitor_class": inhibitor_class,
                "modality_hint": modality_hint,
            }

        tasks = [_process_one(p) for p in raw_proteins]
        results = await asyncio.gather(*tasks, return_exceptions=True)
        for r in results:
            if isinstance(r, dict):
                scored.append(r)

        # Sort by target_score descending, take top 10
        scored.sort(key=lambda x: x["target_score"], reverse=True)
        top10 = scored[:10]

        # Step 4: Assign ranks
        candidates = []
        for i, item in enumerate(top10, start=1):
            candidates.append(TargetCandidate(rank=i, **item))

        # Step 5: Gemini synthesis
        await _step("Generating therapeutic strategy with Gemini...")
        import json
        candidates_json = json.dumps(
            [
                {
                    "rank": c.rank,
                    "gene": c.gene_name,
                    "protein": c.protein_name,
                    "score": c.target_score,
                    "essentiality": c.essentiality_note,
                    "human_identity_pct": c.human_identity_pct,
                    "selectivity_risk": c.selectivity_risk,
                    "known_target": c.known_target,
                    "known_inhibitor_class": c.known_inhibitor_class,
                }
                for c in candidates
            ],
            indent=2,
        )
        from core.gemini_interpreter import synthesize_target_discovery
        synthesis = await synthesize_target_discovery(disease_name, organism_scientific, candidates_json)
        ranked_rationale = synthesis.get("ranked_rationale", "")
        overall_strategy = synthesis.get("overall_strategy", "")

        # Annotate candidates with per-candidate rationale (simple extraction)
        for c in candidates:
            c.gemini_rationale = ranked_rationale  # Full rationale; frontend can display per-card or collapsible

        await _step("Target discovery analysis complete.")
        return TargetDiscoveryReport(
            query=query,
            disease_name=disease_name,
            organism_scientific=organism_scientific,
            organism_common=organism_common,
            candidates=candidates,
            context_summary=context_summary,
            gemini_strategy=overall_strategy,
        )

    except Exception as exc:
        _log.error("run_target_discovery failed: %s", exc, exc_info=True)
        return TargetDiscoveryReport(
            query=query,
            disease_name="unknown",
            organism_scientific="unknown",
            organism_common="unknown",
            context_summary="",
            gemini_strategy="Target discovery analysis failed. Please try again.",
        )
