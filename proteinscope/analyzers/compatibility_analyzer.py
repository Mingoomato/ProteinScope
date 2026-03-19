"""Cross-species functional substitution / de novo expression analyzer.

Two analysis modes:
1. Ortholog substitution — "Can I use penguin POLD1 instead of human POLD1 in HeLa cells?"
   Requires both organisms to have the gene; compares protein sequences, active sites, structure.

2. De novo expression — "Can I express plant RBCL in human epithelial cells?"
   Target organism has NO ortholog; analyzes whether the source gene can be expressed
   heterologously: CDS codon bias, GC content, subcellular targeting signals, organelle
   requirements, post-translational modification compatibility.

Workflow:
1.  Resolve source protein from UniProt
2.  Probe target organism for ortholog (strict — no fallback)
3a. [Ortholog mode] Protein sequence alignment (BLOSUM62) + active site conservation
3b. [De novo mode]  Extract source CDS from NCBI; compute GC content
4.  AlphaFold TM-score (ortholog mode only)
5.  Fetch source CDS from NCBI → Codon Adaptation Index in target organism
6.  Extract subcellular localization + organelle-targeting signals from UniProt
7.  Temperature optima comparison (if available)
8.  Feed ALL computed facts to Gemini 2.5 Pro for synthesis, with mode-aware prompting
"""

from __future__ import annotations

import asyncio
import json as _json
from typing import Optional

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Report model
# ---------------------------------------------------------------------------

class CompatibilityReport(BaseModel):
    source_gene: str
    source_organism: str
    target_gene: str
    target_organism: str
    target_cell: Optional[str] = None

    # Analysis mode
    is_de_novo_expression: bool = False  # True when target has no ortholog

    # Protein-level metrics (ortholog mode)
    sequence_identity_pct: Optional[float] = None
    sequence_similarity_pct: Optional[float] = None
    active_site_conservation_pct: Optional[float] = None
    tm_score: Optional[float] = None
    temperature_delta_celsius: Optional[float] = None
    alignment_details: Optional[str] = None

    # CDS-level metrics (relevant in both modes, critical in de novo mode)
    codon_adaptation_index: Optional[float] = None
    source_gc_content: Optional[float] = None       # GC% of source CDS
    source_cds_length: Optional[int] = None         # CDS length in bp

    # Subcellular biology (critical for de novo feasibility)
    source_subcellular_locations: list[str] = Field(default_factory=list)
    source_targeting_signals: list[str] = Field(default_factory=list)  # transit/signal peptides
    source_organelle: Optional[str] = None           # if organelle-specific (chloroplast, mitochondria)

    # Verdict from Gemini
    verdict: str = ""
    confidence: str = ""
    reasoning: str = ""
    caveats: list[str] = Field(default_factory=list)
    recommendations: list[str] = Field(default_factory=list)

    # UniProt accessions
    source_uniprot_id: Optional[str] = None
    target_uniprot_id: Optional[str] = None

    # Canonical protein sequences (for frontend alignment visualizer)
    source_sequence: Optional[str] = None
    target_sequence: Optional[str] = None

    # Foundational references
    references: list[dict] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# Curated foundational references for cross-species compatibility analysis
# ---------------------------------------------------------------------------

_COMPAT_REFERENCES: list[dict] = [
    {
        "pubmed_id": "2231712",
        "doi": "10.1016/S0022-2836(05)80360-2",
        "title": "Basic local alignment search tool",
        "authors": ["Altschul SF", "Gish W", "Miller W", "Myers EW", "Lipman DJ"],
        "journal": "Journal of Molecular Biology",
        "year": 1990,
        "apa_citation": "Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403–410. https://doi.org/10.1016/S0022-2836(05)80360-2",
    },
    {
        "pubmed_id": "36927031",
        "doi": "10.1126/science.ade2574",
        "title": "Evolutionary-scale prediction of atomic-level protein structure with a language model",
        "authors": ["Lin Z", "Akin H", "Rao R", "Hie B", "Zhu Z", "Lu W", "Rives A"],
        "journal": "Science",
        "year": 2023,
        "apa_citation": "Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., & Rives, A. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. Science, 379(6637), 1123–1130. https://doi.org/10.1126/science.ade2574",
    },
    {
        "pubmed_id": "36408920",
        "doi": "10.1093/nar/gkac1052",
        "title": "UniProt: the Universal Protein Knowledgebase in 2023",
        "authors": ["UniProt Consortium"],
        "journal": "Nucleic Acids Research",
        "year": 2023,
        "apa_citation": "UniProt Consortium. (2023). UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Research, 51(D1), D523–D531. https://doi.org/10.1093/nar/gkac1052",
    },
    {
        "pubmed_id": "3547335",
        "doi": "10.1093/nar/15.3.1281",
        "title": "The codon Adaptation Index — a measure of directional synonymous codon usage bias, and its potential applications",
        "authors": ["Sharp PM", "Li WH"],
        "journal": "Nucleic Acids Research",
        "year": 1987,
        "apa_citation": "Sharp, P. M., & Li, W. H. (1987). The codon Adaptation Index — a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research, 15(3), 1281–1295. https://doi.org/10.1093/nar/15.3.1281",
    },
    {
        "pubmed_id": "25950237",
        "doi": "10.1093/nar/gkv315",
        "title": "KEGG as a reference resource for gene and protein annotation",
        "authors": ["Kanehisa M", "Sato Y", "Kawashima M", "Furumichi M", "Tanabe M"],
        "journal": "Nucleic Acids Research",
        "year": 2016,
        "apa_citation": "Kanehisa, M., Sato, Y., Kawashima, M., Furumichi, M., & Tanabe, M. (2016). KEGG as a reference resource for gene and protein annotation. Nucleic Acids Research, 44(D1), D457–D462. https://doi.org/10.1093/nar/gkv315",
    },
]


# ---------------------------------------------------------------------------
# Helpers — sequence extraction
# ---------------------------------------------------------------------------

def _extract_sequence(entry: dict) -> str:
    return entry.get("sequence", {}).get("value", "")


def _active_site_dicts(entry: dict) -> list[dict]:
    sites = []
    for feat in entry.get("features", []):
        if feat.get("type", "").lower() == "active site":
            loc = feat.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            if start is not None and end is not None:
                sites.append({
                    "start": int(start) - 1,
                    "end": int(end),
                    "description": feat.get("description", "active site"),
                })
    return sites


def _extract_temperature(entry: dict) -> Optional[float]:
    import re
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "BIOPHYSICOCHEMICAL PROPERTIES":
            for t in comment.get("temperatureDependence", {}).get("texts", []):
                m = re.search(r"optimum.*?(\d+(?:\.\d+)?)\s*[°℃C]", t.get("value", ""), re.IGNORECASE)
                if m:
                    return float(m.group(1))
    return None


# ---------------------------------------------------------------------------
# Helpers — subcellular biology
# ---------------------------------------------------------------------------

def _extract_subcellular_locations(entry: dict) -> list[str]:
    """Return list of subcellular localisation strings from UniProt."""
    locs = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for sub in comment.get("subcellularLocations", []):
                val = sub.get("location", {}).get("value", "")
                if val:
                    locs.append(val)
    return locs


def _extract_targeting_signals(entry: dict) -> tuple[list[str], Optional[str]]:
    """Return (list of signal descriptions, organelle name or None).

    Scans UniProt 'features' for transit peptides, signal peptides, propeptides.
    Also infers organelle (Chloroplast / Mitochondrion / ER / Nucleus …)
    from the subcellular location strings.
    """
    signals: list[str] = []
    organelle: Optional[str] = None

    signal_types = {"transit peptide", "signal peptide", "propeptide"}
    for feat in entry.get("features", []):
        ftype = feat.get("type", "").lower()
        if ftype in signal_types:
            desc = feat.get("description", ftype).strip() or ftype
            signals.append(desc)

    # Infer organelle from subcellular location text
    locations = _extract_subcellular_locations(entry)
    organelle_keywords = {
        "Chloroplast": "Chloroplast",
        "Mitochondri": "Mitochondrion",
        "Nucleus":     "Nucleus",
        "Endoplasmic reticulum": "Endoplasmic reticulum",
        "Golgi":       "Golgi apparatus",
        "Peroxisome":  "Peroxisome",
    }
    for loc in locations:
        for kw, name in organelle_keywords.items():
            if kw.lower() in loc.lower():
                organelle = name
                break
        if organelle:
            break

    return signals, organelle


# ---------------------------------------------------------------------------
# Helpers — CDS / sequence analytics
# ---------------------------------------------------------------------------

def _compute_gc_content(cds: str) -> Optional[float]:
    """GC% of a nucleotide sequence."""
    if not cds:
        return None
    upper = cds.upper()
    total = sum(1 for c in upper if c in "ACGT")
    gc = sum(1 for c in upper if c in "GC")
    return round(100.0 * gc / total, 1) if total > 0 else None


def _compute_cai(cds: str, target_organism: str) -> Optional[float]:
    """Codon Adaptation Index of source CDS in target organism codon table."""
    if not cds or len(cds) < 9:
        return None
    try:
        from Bio.Seq import Seq
        from Bio.Data.CodonTable import standard_dna_table
        import math

        seq = Seq(cds[:len(cds) - (len(cds) % 3)])
        codons: dict[str, int] = {}
        for i in range(0, len(seq) - 2, 3):
            codon = str(seq[i:i+3])
            if len(codon) == 3 and "N" not in codon:
                codons[codon] = codons.get(codon, 0) + 1
        if not codons:
            return None

        aa_to_codons: dict[str, list[str]] = {}
        for codon, aa in standard_dna_table.forward_table.items():
            aa_to_codons.setdefault(aa, []).append(codon)

        w: dict[str, float] = {}
        for aa, syn_codons in aa_to_codons.items():
            counts = [codons.get(c, 0) for c in syn_codons]
            max_count = max(counts) if counts else 0
            for c, cnt in zip(syn_codons, counts):
                w[c] = cnt / max_count if max_count > 0 else 0.0

        log_sum = 0.0
        n = 0
        for i in range(0, len(seq) - 2, 3):
            codon = str(seq[i:i+3])
            wi = w.get(codon)
            if wi and wi > 0:
                log_sum += math.log(wi)
                n += 1
        if n == 0:
            return None
        return round(math.exp(log_sum / n), 4)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# UniProt lookup
# ---------------------------------------------------------------------------

async def _fetch_uniprot_entry(gene: str, organism: str, strict: bool = False) -> Optional[dict]:
    """Search UniProt for a gene in an organism.

    strict=True → only search within the given organism (no any-organism fallback).
    Use strict=True for the *target* to avoid the self-comparison bug where a missing
    ortholog silently returns the source protein.
    """
    try:
        from fetchers import uniprot as uni
        candidates = [
            f"reviewed:true AND gene_exact:{gene}",
            f"reviewed:true AND protein_name:{gene}",
            gene,
        ]
        org_attempts = [organism] if strict else [organism, None]
        for org_try in org_attempts:
            for q in candidates:
                results = await uni.search_protein(q, org_try, size=1)
                if results:
                    accession = results[0]["primaryAccession"]
                    entry = await uni.fetch_by_accession(accession)
                    return entry
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Structure comparison
# ---------------------------------------------------------------------------

async def _compute_tm_score(source_accession: str, target_accession: str) -> Optional[float]:
    try:
        from fetchers.alphafold import download_pdb
        import httpx, numpy as np
        from pathlib import Path

        async with httpx.AsyncClient(timeout=60) as client:
            src_path = await download_pdb(source_accession, client)
            tgt_path = await download_pdb(target_accession, client)

        try:
            import biotite.structure.io.pdb as pdbio

            def _read_ca(path: Path):
                f = pdbio.PDBFile.read(str(path))
                arr = f.get_structure(model=1, extra_fields=[])
                return arr[arr.atom_name == "CA"].coord

            c1 = _read_ca(src_path)
            c2 = _read_ca(tgt_path)
            min_len = min(len(c1), len(c2))
            if min_len == 0:
                return None
            c1, c2 = c1[:min_len], c2[:min_len]
            d0 = max(1.24 * (min_len - 15) ** (1 / 3) - 1.8, 0.5) if min_len > 15 else 0.5
            return round(float(np.mean(1 / (1 + (np.sqrt(np.sum((c1 - c2) ** 2, axis=1)) / d0) ** 2))), 4)
        except Exception:
            pass
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Gemini synthesis — mode-aware
# ---------------------------------------------------------------------------

async def _synthesize_with_gemini(report: CompatibilityReport) -> CompatibilityReport:
    try:
        from core.gemini_interpreter import _call

        facts_lines: list[str] = []

        if report.is_de_novo_expression:
            facts_lines.append(
                f"ANALYSIS MODE: De novo expression "
                f"(no {report.source_gene} ortholog exists in {report.target_organism}; "
                f"the goal is to EXPRESS the {report.source_organism} gene heterologously)"
            )
        else:
            facts_lines.append(
                f"ANALYSIS MODE: Ortholog substitution "
                f"({report.target_organism} has a {report.source_gene} ortholog; "
                f"evaluating functional interchangeability)"
            )

        # Protein-level evidence (ortholog mode)
        if report.sequence_identity_pct is not None:
            facts_lines.append(f"- Protein sequence identity: {report.sequence_identity_pct:.1f}%")
        if report.sequence_similarity_pct is not None:
            facts_lines.append(f"- Protein sequence similarity: {report.sequence_similarity_pct:.1f}%")
        if report.active_site_conservation_pct is not None:
            facts_lines.append(f"- Active site conservation: {report.active_site_conservation_pct:.1f}%")
        if report.tm_score is not None:
            facts_lines.append(f"- AlphaFold structural TM-score: {report.tm_score:.3f}")
        if report.temperature_delta_celsius is not None:
            facts_lines.append(f"- Temperature optimum delta: {report.temperature_delta_celsius:+.1f} °C")

        # CDS / codon evidence
        if report.codon_adaptation_index is not None:
            cai_interp = (
                "good codon match" if report.codon_adaptation_index > 0.7
                else "moderate codon mismatch" if report.codon_adaptation_index > 0.4
                else "severe codon mismatch — low expression expected"
            )
            facts_lines.append(
                f"- Codon Adaptation Index in {report.target_organism}: "
                f"{report.codon_adaptation_index:.3f} ({cai_interp})"
            )
        if report.source_gc_content is not None:
            facts_lines.append(f"- Source CDS GC content: {report.source_gc_content:.1f}%")
            # Human genome avg GC ~41%; plant genomes often differ
            human_gc = 41.0
            delta_gc = abs(report.source_gc_content - human_gc)
            if delta_gc > 15:
                facts_lines.append(
                    f"  → GC content diverges {delta_gc:.1f}% from human average (~41%) — "
                    f"may affect RNA stability and translation efficiency"
                )
        if report.source_cds_length is not None:
            facts_lines.append(f"- Source CDS length: {report.source_cds_length} bp")

        # Subcellular biology
        if report.source_subcellular_locations:
            facts_lines.append(
                f"- Native subcellular localization: {', '.join(report.source_subcellular_locations)}"
            )
        if report.source_organelle:
            facts_lines.append(
                f"- Organelle requirement: {report.source_organelle} "
                f"({'ABSENT in ' + report.target_organism if report.source_organelle == 'Chloroplast' else 'present in ' + report.target_organism})"
            )
        if report.source_targeting_signals:
            facts_lines.append(
                f"- Targeting/signal sequences detected: {', '.join(report.source_targeting_signals)}"
            )

        if not facts_lines:
            facts_lines.append("- No quantitative metrics were computed (sequence databases unavailable).")

        facts_block = "\n".join(facts_lines)
        cell_ctx = f"in {report.target_cell} cells" if report.target_cell else f"in {report.target_organism} cells"

        if report.is_de_novo_expression:
            question = (
                f"Can {report.source_organism} {report.source_gene} be expressed de novo {cell_ctx} "
                f"and perform its biological function? "
                f"(Note: {report.target_organism} has NO endogenous {report.source_gene} ortholog.)"
            )
            mode_guidance = (
                "IMPORTANT: This is a DE NOVO EXPRESSION analysis, not ortholog substitution.\n"
                "Address:\n"
                "  1. Fundamental biological barriers (missing organelles, absent cofactors, wrong cellular compartment)\n"
                "  2. Codon optimisation requirements based on CAI and GC content\n"
                "  3. What modifications to the gene/construct would be needed\n"
                "  4. Whether the protein's function is even relevant in target cells\n"
                "  5. Historical precedent (has this ever been attempted / succeeded?)\n"
                "Be biologically rigorous — do not say 'compatible' just because sequence identity looks high.\n"
                "If a critical organelle (e.g. chloroplast) is absent, that alone makes function impossible.\n"
            )
        else:
            question = (
                f"Can {report.source_organism} {report.source_gene} functionally substitute "
                f"for {report.target_organism} {report.source_gene} {cell_ctx}?"
            )
            mode_guidance = (
                "Thresholds for guidance:\n"
                "  sequence_identity >= 90% → likely functional conservation\n"
                "  sequence_identity 70–89% → possible, verify active sites carefully\n"
                "  sequence_identity < 70%  → high risk of incompatibility\n"
                "  TM-score >= 0.9          → highly similar structure\n"
                "  active_site_conservation < 100% → flag specific residue changes\n"
                "  CAI < 0.4               → codon optimisation strongly recommended\n"
            )

        prompt = (
            "You are a molecular and cellular biology expert.\n\n"
            f"Question: {question}\n\n"
            "Computed evidence:\n"
            f"{facts_block}\n\n"
            f"{mode_guidance}\n"
            "Return a raw JSON object with these fields:\n"
            "  verdict: one of \"Compatible\", \"Likely compatible\", \"Conditional\", "
            "\"Incompatible\", or \"Unknown\"\n"
            "  confidence: \"High\", \"Moderate\", or \"Low\"\n"
            "  reasoning: 3-5 sentences grounded strictly in the computed evidence above — "
            "explain the biology of WHY, not just the numbers\n"
            "  caveats: JSON array of specific biological concerns (minimum 3)\n"
            "  recommendations: JSON array of next experimental steps (minimum 3)\n\n"
            "Return ONLY the raw JSON, no markdown fences:"
        )

        raw = await _call(prompt)
        if not raw:
            report.verdict = "Unknown"
            report.confidence = "Low"
            report.reasoning = "AI synthesis unavailable — GEMINI_API_KEY not set or quota exceeded."
            return report

        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)

        valid_verdicts = {"Compatible", "Likely compatible", "Conditional", "Incompatible", "Unknown"}
        verdict = str(data.get("verdict", "Unknown")).strip()
        report.verdict = verdict if verdict in valid_verdicts else "Unknown"

        valid_conf = {"High", "Moderate", "Low"}
        conf = str(data.get("confidence", "Low")).strip()
        report.confidence = conf if conf in valid_conf else "Low"

        report.reasoning = str(data.get("reasoning", "")).strip()

        caveats = data.get("caveats", [])
        report.caveats = [str(c).strip() for c in caveats if c] if isinstance(caveats, list) else []

        recs = data.get("recommendations", [])
        report.recommendations = [str(r).strip() for r in recs if r] if isinstance(recs, list) else []

    except Exception:
        if not report.verdict:
            report.verdict = "Unknown"
            report.confidence = "Low"
            report.reasoning = "AI synthesis failed due to an unexpected error."

    return report


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_compatibility_analysis(
    source_gene: str,
    source_organism: str,
    target_organism: str,
    target_cell: Optional[str] = None,
    step_cb=None,
) -> CompatibilityReport:
    """Run full cross-species compatibility / de novo expression analysis.

    Automatically detects whether an ortholog exists in the target organism and
    switches between ortholog-substitution and de-novo-expression analysis modes.
    """
    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    report = CompatibilityReport(
        source_gene=source_gene.upper(),
        source_organism=source_organism,
        target_gene=source_gene.upper(),
        target_organism=target_organism,
        target_cell=target_cell,
    )

    # ── 1. Resolve UniProt entries ──────────────────────────────────────────
    await _step(f"[1/7] Looking up {source_gene} in UniProt ({source_organism} & {target_organism})…")
    source_entry, target_entry = await asyncio.gather(
        _fetch_uniprot_entry(source_gene, source_organism, strict=False),
        _fetch_uniprot_entry(source_gene, target_organism, strict=True),
        return_exceptions=True,
    )
    if isinstance(source_entry, Exception):
        source_entry = None
    if isinstance(target_entry, Exception):
        target_entry = None

    if isinstance(source_entry, dict):
        report.source_uniprot_id = source_entry.get("primaryAccession") or source_entry.get("uniProtkbId")
    if isinstance(target_entry, dict):
        report.target_uniprot_id = target_entry.get("primaryAccession") or target_entry.get("uniProtkbId")

    # Self-comparison guard
    if (report.source_uniprot_id and report.target_uniprot_id
            and report.source_uniprot_id == report.target_uniprot_id):
        report.target_uniprot_id = None
        target_entry = None

    # Determine analysis mode
    report.is_de_novo_expression = (
        isinstance(source_entry, dict) and not isinstance(target_entry, dict)
    )

    source_seq = _extract_sequence(source_entry) if isinstance(source_entry, dict) else ""
    target_seq = _extract_sequence(target_entry) if isinstance(target_entry, dict) else ""

    # Store sequences for frontend alignment visualizer (cap at 2000 aa to avoid bloat)
    if source_seq:
        report.source_sequence = source_seq[:2000]
    if target_seq:
        report.target_sequence = target_seq[:2000]

    # ── 2. Subcellular biology (always — critical for de novo mode) ─────────
    await _step("[2/7] Extracting subcellular localization & targeting signals…")
    if isinstance(source_entry, dict):
        report.source_subcellular_locations = _extract_subcellular_locations(source_entry)
        signals, organelle = _extract_targeting_signals(source_entry)
        report.source_targeting_signals = signals
        report.source_organelle = organelle

    # ── 3. Protein sequence alignment (ortholog mode only) ──────────────────
    await _step("[3/7] Running BLOSUM62 sequence alignment…")
    if source_seq and target_seq:
        try:
            from analyzers.sequence_aligner import align_pairwise, align_at_annotated_sites
            aln = align_pairwise(source_seq, target_seq)
            if aln:
                report.sequence_identity_pct = aln.identity_pct
                report.sequence_similarity_pct = aln.similarity_pct
                report.alignment_details = (
                    f"Identity: {aln.identity_pct:.1f}%, "
                    f"Similarity: {aln.similarity_pct:.1f}%, "
                    f"Gaps: {aln.gap_pct:.1f}%"
                )
            active_sites_src = _active_site_dicts(source_entry) if isinstance(source_entry, dict) else []
            if active_sites_src:
                site_results = align_at_annotated_sites(source_seq, target_seq, active_sites_src)
                if site_results:
                    conserved = sum(1 for s in site_results if s.get("conserved"))
                    total = len(site_results)
                    report.active_site_conservation_pct = round(100.0 * conserved / total, 1) if total else None
        except Exception:
            pass
    elif report.is_de_novo_expression:
        await _step("[3/7] De novo mode — skipping ortholog alignment (no target sequence)")

    # ── 4. AlphaFold structural comparison (ortholog mode only) ────────────
    await _step("[4/7] Downloading AlphaFold structures & computing TM-score…")
    if report.source_uniprot_id and report.target_uniprot_id:
        try:
            report.tm_score = await _compute_tm_score(report.source_uniprot_id, report.target_uniprot_id)
        except Exception:
            pass

    # ── 5. Fetch source CDS & compute GC content + CAI ─────────────────────
    await _step("[5/7] Fetching source CDS from NCBI → GC content & Codon Adaptation Index…")
    cds_seq: Optional[str] = None
    if isinstance(source_entry, dict):
        try:
            from fetchers import uniprot as uni
            ncbi_gene_id = uni.get_ncbi_gene_id(source_entry)
            if ncbi_gene_id:
                from fetchers.ncbi_gene import fetch_cds_for_gene
                cds_result = await fetch_cds_for_gene(ncbi_gene_id)
                if isinstance(cds_result, tuple):
                    cds_seq = cds_result[0]
                elif isinstance(cds_result, str):
                    cds_seq = cds_result
        except Exception:
            pass

    if cds_seq:
        report.source_gc_content = _compute_gc_content(cds_seq)
        report.source_cds_length = len(cds_seq)
        report.codon_adaptation_index = _compute_cai(cds_seq, target_organism)
    else:
        # Fall back to protein sequence as proxy for CAI (less accurate)
        report.codon_adaptation_index = _compute_cai(source_seq, target_organism)

    # ── 6. Temperature optima comparison ───────────────────────────────────
    await _step("[6/7] Comparing temperature optima (BRENDA)…")
    try:
        from fetchers.brenda import get_ec_number, fetch_brenda_temperature, extract_temperature_from_uniprot
        src_temp = tgt_temp = None
        if isinstance(source_entry, dict):
            ec = get_ec_number(source_entry)
            if ec:
                src_temp = await fetch_brenda_temperature(ec, source_organism)
            if src_temp is None:
                src_temp = extract_temperature_from_uniprot(source_entry)
        if isinstance(target_entry, dict):
            ec = get_ec_number(target_entry)
            if ec:
                tgt_temp = await fetch_brenda_temperature(ec, target_organism)
            if tgt_temp is None:
                tgt_temp = extract_temperature_from_uniprot(target_entry)
        if src_temp is not None and tgt_temp is not None:
            report.temperature_delta_celsius = round(src_temp - tgt_temp, 1)
    except Exception:
        pass

    # ── 7. Gemini synthesis ────────────────────────────────────────────────
    await _step("[7/7] Gemini 2.5 Pro synthesizing verdict…")
    report = await _synthesize_with_gemini(report)

    report.references = _COMPAT_REFERENCES
    return report


async def run_compatibility_with_impact(
    source_gene: str,
    source_organism: str,
    target_organism: str,
    target_cell: Optional[str] = None,
    step_cb=None,
) -> tuple["CompatibilityReport", "InsertionImpactReport"]:
    """Run compatibility analysis AND metabolic insertion impact analysis together.

    Returns (CompatibilityReport, InsertionImpactReport).
    The insertion impact analysis reuses data already fetched by compatibility analysis.
    """
    from analyzers.insertion_impact_analyzer import (
        InsertionImpactReport,
        run_insertion_impact_analysis,
    )

    compat = await run_compatibility_analysis(
        source_gene=source_gene,
        source_organism=source_organism,
        target_organism=target_organism,
        target_cell=target_cell,
        step_cb=step_cb,
    )

    if step_cb:
        try:
            await step_cb("Running metabolic pathway conflict analysis…")
        except Exception:
            pass

    impact = await run_insertion_impact_analysis(
        source_gene=source_gene,
        source_organism=source_organism,
        target_organism=target_organism,
        source_uniprot_id=compat.source_uniprot_id,
        source_ec_number=None,  # resolved internally from UniProt
        source_organelle=compat.source_organelle,
        step_cb=step_cb,
    )

    return compat, impact
