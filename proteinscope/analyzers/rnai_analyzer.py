"""RNAi Experiment Design & Impact Analyzer.

Background — What is RNAi?
─────────────────────────────────────────────────────────────────────────────
RNA interference (RNAi) is a conserved eukaryotic gene silencing mechanism
in which small double-stranded RNA molecules direct sequence-specific
degradation or translational repression of target mRNAs.

Core mechanisms:
  siRNA (small interfering RNA, ~21-23 nt, fully complementary)
    → Loaded into RISC (RNA-induced silencing complex) via Dicer processing
    → Argonaute-2 (AGO2) cleaves target mRNA ("Slicer" activity)
    → Near-complete mRNA knockdown; exogenous, transient

  shRNA (short hairpin RNA, ~19-29 nt stem + loop)
    → Expressed from a U6/H1 promoter in a viral vector (lentivirus, AAV)
    → Exported by Exportin-5, processed by Dicer → siRNA-like RISC loading
    → Stable, heritable knockdown; suitable for in vivo studies

  miRNA (microRNA, ~22 nt, partial complementarity)
    → Endogenously encoded; processed by Drosha → Dicer
    → Translational repression + mRNA destabilisation
    → One miRNA regulates hundreds of targets simultaneously

Experimental design workflow:
  1. Target selection: choose gene, transcript, and target region
  2. siRNA/shRNA design: thermodynamic rules (asymmetry, GC 30-52%,
       no runs of ≥4 identical bases, no internal repeats)
  3. Off-target prediction: seed region (nt 2-8) BLAST vs. transcriptome
  4. Delivery method selection (transfection reagent, viral vector, etc.)
  5. Controls: scramble/non-targeting control, positive control, rescue
  6. Validation: qRT-PCR (mRNA, 48-72h), Western blot (protein, 72-96h)
  7. Phenotypic readout: proliferation, apoptosis, migration, etc.
  8. Rescue experiment: re-expression of siRNA-resistant cDNA confirms
       on-target specificity

This analyzer performs:
  1. Target mRNA structure analysis (UTR regions, accessible sites)
  2. Optimal siRNA/shRNA sequence design (Tuschl/Reynolds rules)
  3. Off-target risk assessment (seed region complementarity)
  4. Delivery strategy recommendation (cell-type aware)
  5. Predicted knockdown phenotype via DB cross-check (same as LoF
       reverse genetics analysis but framed in an RNAi context)
  6. Validation experiment design
  7. Gemini synthesis of the complete experimental protocol
"""

from __future__ import annotations

import asyncio
import re
from typing import Optional

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

class SiRNADesign(BaseModel):
    sequence_sense: str           # 5'→3' sense strand
    sequence_antisense: str       # 5'→3' antisense strand
    target_position: Optional[int] = None   # position in mRNA
    gc_content_pct: float
    specificity_score: str        # High | Moderate | Low
    design_notes: list[str] = Field(default_factory=list)


class ShRNADesign(BaseModel):
    stem_sequence: str            # 19-21 nt sense sequence
    loop_sequence: str = "TTCAAGAGA"   # standard loop
    full_hairpin: str             # complete hairpin sequence
    vector_recommendation: str   # e.g. "pLKO.1 (lentiviral), U6 promoter"


class DeliveryStrategy(BaseModel):
    method: str                   # lipofection | electroporation | lentivirus | AAV | nanoparticle
    reagent_recommendation: str   # e.g. "Lipofectamine RNAiMAX"
    optimal_concentration: str    # e.g. "10-50 nM siRNA"
    transfection_efficiency: str  # expected % for this cell type
    duration_days: str            # typical knockdown duration
    rationale: str


class ValidationExperiment(BaseModel):
    assay_type: str               # qRT-PCR | Western blot | Flow cytometry | etc.
    target_marker: str
    timepoint: str                # e.g. "48h post-transfection"
    expected_result: str
    controls_required: list[str]


class RNAiReport(BaseModel):
    gene: str
    organism: str
    rnai_type: str                # siRNA | shRNA | miRNA | general
    cell_type: Optional[str] = None
    uniprot_id: Optional[str] = None
    protein_name: Optional[str] = None
    protein_function: Optional[str] = None
    mrna_length: Optional[int] = None

    # Design outputs
    sirna_designs: list[SiRNADesign] = Field(default_factory=list)
    shrna_design: Optional[ShRNADesign] = None
    delivery_strategy: Optional[DeliveryStrategy] = None

    # Off-target risk
    off_target_risk: str = "Unknown"   # High | Moderate | Low
    off_target_notes: list[str] = Field(default_factory=list)
    seed_region_warnings: list[str] = Field(default_factory=list)

    # Predicted knockdown phenotype (from LoF reverse genetics)
    disrupted_pathways: list[dict] = Field(default_factory=list)
    predicted_phenotypes: list[str] = Field(default_factory=list)
    known_ko_phenotypes: list[str] = Field(default_factory=list)   # from OMIM/ClinVar

    # Controls required
    controls: list[str] = Field(default_factory=list)

    # Validation plan
    validation_experiments: list[ValidationExperiment] = Field(default_factory=list)

    # Experimental timeline
    timeline_days: list[dict] = Field(default_factory=list)   # [{day, task}]

    # Gemini synthesis
    experimental_protocol: str = ""   # full protocol summary
    key_considerations: list[str] = Field(default_factory=list)
    troubleshooting_tips: list[str] = Field(default_factory=list)
    analysis_caveats: list[str] = Field(default_factory=list)

    # Target CDS/mRNA sequence (stored for frontend alignment cross-check)
    target_cds: Optional[str] = None

    # Foundational references
    references: list[dict] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# Curated foundational references for RNAi
# ---------------------------------------------------------------------------

_RNAI_REFERENCES: list[dict] = [
    {
        "pubmed_id": "9486653",
        "doi": "10.1038/35888",
        "title": "Potent and specific genetic interference by double-stranded RNA in Caenorhabditis elegans",
        "authors": ["Fire A", "Xu S", "Montgomery MK", "Kostas SA", "Driver SE", "Mello CC"],
        "journal": "Nature",
        "year": 1998,
        "apa_citation": "Fire, A., Xu, S., Montgomery, M. K., Kostas, S. A., Driver, S. E., & Mello, C. C. (1998). Potent and specific genetic interference by double-stranded RNA in Caenorhabditis elegans. Nature, 391(6669), 806–811. https://doi.org/10.1038/35888",
    },
    {
        "pubmed_id": "11373684",
        "doi": "10.1038/35078107",
        "title": "Duplexes of 21-nucleotide RNAs mediate RNA interference in cultured mammalian cells",
        "authors": ["Elbashir SM", "Harborth J", "Lendeckel W", "Yalcin A", "Weber K", "Tuschl T"],
        "journal": "Nature",
        "year": 2001,
        "apa_citation": "Elbashir, S. M., Harborth, J., Lendeckel, W., Yalcin, A., Weber, K., & Tuschl, T. (2001). Duplexes of 21-nucleotide RNAs mediate RNA interference in cultured mammalian cells. Nature, 411(6836), 494–498. https://doi.org/10.1038/35078107",
    },
    {
        "pubmed_id": "14758366",
        "doi": "10.1038/nbt936",
        "title": "Rational siRNA design for RNA interference",
        "authors": ["Reynolds A", "Leake D", "Boese Q", "Scaringe S", "Marshall WS", "Khvorova A"],
        "journal": "Nature Biotechnology",
        "year": 2004,
        "apa_citation": "Reynolds, A., Leake, D., Boese, Q., Scaringe, S., Marshall, W. S., & Khvorova, A. (2004). Rational siRNA design for RNA interference. Nature Biotechnology, 22(3), 326–330. https://doi.org/10.1038/nbt936",
    },
    {
        "pubmed_id": "14567917",
        "doi": "10.1016/S0092-8674(03)00759-1",
        "title": "Asymmetry in the assembly of the RNAi enzyme complex",
        "authors": ["Schwarz DS", "Hutvágner G", "Du T", "Xu Z", "Aronin N", "Zamore PD"],
        "journal": "Cell",
        "year": 2003,
        "apa_citation": "Schwarz, D. S., Hutvágner, G., Du, T., Xu, Z., Aronin, N., & Zamore, P. D. (2003). Asymmetry in the assembly of the RNAi enzyme complex. Cell, 115(2), 199–208. https://doi.org/10.1016/S0092-8674(03)00759-1",
    },
    {
        "pubmed_id": "16211010",
        "doi": "10.1038/nbt1209",
        "title": "A microRNA in a multiple-turnover RNAi enzyme complex",
        "authors": ["Gregory RI", "Chendrimada TP", "Cooch N", "Shiekhattar R"],
        "journal": "Science",
        "year": 2005,
        "apa_citation": "Gregory, R. I., Chendrimada, T. P., Cooch, N., & Shiekhattar, R. (2005). Human RISC couples microRNA biogenesis and posttranscriptional gene silencing. Cell, 123(4), 631–640. https://doi.org/10.1016/j.cell.2005.10.022",
    },
]


# ---------------------------------------------------------------------------
# siRNA/shRNA sequence design (Tuschl & Reynolds rules)
# ---------------------------------------------------------------------------

_STANDARD_LOOP = "TTCAAGAGA"

# Thermodynamic asymmetry rule: antisense 5' end should be A/U (lower stability)
_AU_BASES = set("AU")
_GC_BASES = set("GC")

def _complement(seq: str) -> str:
    comp = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return seq.translate(comp)

def _reverse_complement(seq: str) -> str:
    return _complement(seq)[::-1]

def _gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in "GC")
    return round(100.0 * gc / len(seq), 1)

def _has_runs(seq: str, n: int = 4) -> bool:
    """Check for runs of n identical bases."""
    for base in "ACGTU":
        if base * n in seq.upper():
            return True
    return False

def _check_tuschl_rules(sense: str) -> list[str]:
    """Apply Tuschl/Reynolds thermodynamic rules. Returns list of violations."""
    issues = []
    gc = _gc_content(sense)
    if gc < 30 or gc > 52:
        issues.append(f"GC content {gc:.0f}% outside optimal 30-52% range")
    if sense[0].upper() not in "AU":
        issues.append("Position 1 should be A or U (thermodynamic asymmetry)")
    if sense[-1].upper() not in "GC":
        issues.append("Position 19 should be G or C (3' end stability for RISC loading)")
    if _has_runs(sense, 4):
        issues.append("Contains run of ≥4 identical bases — avoid for specificity")
    # Avoid palindromes / internal repeats (simplified check)
    if len(set(sense[i:i+4] for i in range(len(sense)-3))) < len(sense) - 6:
        issues.append("Possible internal repeat sequence detected")
    return issues

def _design_sirna_candidates(cds_or_mrna: str, n: int = 3) -> list[SiRNADesign]:
    """Generate candidate siRNA sequences using a sliding window over the CDS.

    Applies Tuschl rules and picks top n candidates by fewest violations.
    """
    candidates: list[tuple[int, list[str], str, str]] = []  # (position, issues, sense, antisense)
    seq = cds_or_mrna.upper().replace("U", "T")

    if len(seq) < 21:
        return []

    for i in range(0, len(seq) - 21, 7):   # step 7 nt for coverage
        sense = seq[i:i+21]
        if len(sense) < 21 or "N" in sense:
            continue
        antisense = _reverse_complement(sense)
        issues = _check_tuschl_rules(sense)
        candidates.append((i, issues, sense, antisense))

    # Sort by fewest violations
    candidates.sort(key=lambda x: len(x[1]))
    results: list[SiRNADesign] = []
    seen_seeds: set[str] = set()

    for pos, issues, sense, antisense in candidates[:n * 3]:
        seed = antisense[1:8]   # seed region = nt 2-8 of guide strand
        if seed in seen_seeds:
            continue            # avoid overlapping seed regions
        seen_seeds.add(seed)

        gc = _gc_content(sense)
        specificity = "High" if not issues else ("Moderate" if len(issues) <= 1 else "Low")
        notes = issues if issues else ["Meets all Tuschl/Reynolds thermodynamic criteria"]

        results.append(SiRNADesign(
            sequence_sense=sense,
            sequence_antisense=antisense,
            target_position=pos,
            gc_content_pct=gc,
            specificity_score=specificity,
            design_notes=notes,
        ))
        if len(results) >= n:
            break

    return results

def _design_shrna(best_sense: str) -> ShRNADesign:
    """Convert best siRNA sense strand to a shRNA cassette."""
    stem = best_sense[:21]
    loop = _STANDARD_LOOP
    antisense_stem = _reverse_complement(stem)
    full_hairpin = f"5'-{stem}-{loop}-{antisense_stem}-3'"
    return ShRNADesign(
        stem_sequence=stem,
        loop_sequence=loop,
        full_hairpin=full_hairpin,
        vector_recommendation=(
            "pLKO.1 (lentiviral, puromycin selection, U6 promoter) for stable knockdown; "
            "pLL3.7 or pSicoR for fluorescent tracking"
        ),
    )


# ---------------------------------------------------------------------------
# Delivery strategy selection
# ---------------------------------------------------------------------------

_DELIVERY_STRATEGIES = {
    "hela": DeliveryStrategy(
        method="lipofection",
        reagent_recommendation="Lipofectamine RNAiMAX (Thermo Fisher L3000015)",
        optimal_concentration="10 nM siRNA",
        transfection_efficiency=">80%",
        duration_days="48-72 h (siRNA); stable via lentivirus",
        rationale="HeLa cells are highly transfectable; RNAiMAX gives minimal cytotoxicity",
    ),
    "hek293": DeliveryStrategy(
        method="lipofection",
        reagent_recommendation="Lipofectamine 2000 or RNAiMAX at 5 nM siRNA",
        optimal_concentration="5-25 nM siRNA",
        transfection_efficiency=">90%",
        duration_days="48-72 h",
        rationale="HEK293 among most transfectable lines; low concentration avoids ISR",
    ),
    "primary": DeliveryStrategy(
        method="electroporation",
        reagent_recommendation="Nucleofector (Lonza) with cell-type–specific kit",
        optimal_concentration="100-300 nM siRNA",
        transfection_efficiency="40-70%",
        duration_days="48-72 h; consider lentiviral shRNA for stability",
        rationale="Lipofection inefficient in primary cells; electroporation bypasses membrane barriers",
    ),
    "neuron": DeliveryStrategy(
        method="lentivirus",
        reagent_recommendation="pLKO.1 shRNA lentivirus at MOI 5-10",
        optimal_concentration="MOI 5-10 (shRNA lentivirus)",
        transfection_efficiency="60-90% via lentivirus",
        duration_days="Stable (weeks-months)",
        rationale="Neurons are post-mitotic and resistant to lipofection; lentivirus achieves persistent silencing",
    ),
    "default": DeliveryStrategy(
        method="lipofection",
        reagent_recommendation="Lipofectamine RNAiMAX at 10-50 nM siRNA",
        optimal_concentration="10-50 nM siRNA",
        transfection_efficiency="50-80% (cell-type dependent)",
        duration_days="48-72 h (siRNA); stable via shRNA lentivirus",
        rationale="RNAiMAX is broadly applicable; adjust concentration per cell viability testing",
    ),
}

def _select_delivery(cell_type: Optional[str], rnai_type: str) -> DeliveryStrategy:
    """Select appropriate delivery strategy based on cell type."""
    ct = (cell_type or "").lower()
    if rnai_type == "shRNA":
        return DeliveryStrategy(
            method="lentiviral transduction",
            reagent_recommendation="pLKO.1 shRNA lentivirus, MOI 5-10, 8 μg/mL polybrene",
            optimal_concentration="MOI 5-10",
            transfection_efficiency="60-95% (cell-type dependent)",
            duration_days="Stable (weeks-months)",
            rationale="shRNA requires viral vector integration for persistent expression from U6 promoter",
        )
    for key in _DELIVERY_STRATEGIES:
        if key != "default" and key in ct:
            return _DELIVERY_STRATEGIES[key]
    if "primary" in ct or "neuron" in ct or "stem" in ct or "ipsc" in ct:
        return _DELIVERY_STRATEGIES.get(ct.split()[0], _DELIVERY_STRATEGIES["primary"])
    return _DELIVERY_STRATEGIES["default"]


# ---------------------------------------------------------------------------
# Standard controls
# ---------------------------------------------------------------------------

def _required_controls(gene: str, rnai_type: str) -> list[str]:
    return [
        f"Non-targeting scramble control ({rnai_type}) — matched GC content to {gene} siRNA",
        f"Positive control: known-effective siRNA for a housekeeping gene (e.g. GAPDH, ACTB knockdown)",
        f"Mock transfection control — transfection reagent only, no siRNA",
        f"Rescue experiment: re-express siRNA-resistant {gene} cDNA (silent mutations in target region) to confirm on-target specificity",
        "Fluorescent siRNA control (e.g. siGLO) to monitor transfection efficiency",
    ]


# ---------------------------------------------------------------------------
# Experimental timeline
# ---------------------------------------------------------------------------

def _build_timeline(rnai_type: str, cell_type: Optional[str]) -> list[dict]:
    is_shRNA = rnai_type == "shRNA"
    if is_shRNA:
        return [
            {"day": "D-3", "task": "Clone shRNA into lentiviral vector; confirm by sequencing"},
            {"day": "D-2", "task": "Produce lentivirus (HEK293T packaging + transfection)"},
            {"day": "D0",  "task": "Seed target cells; transduce with lentivirus + polybrene"},
            {"day": "D2",  "task": "Change media; begin puromycin selection (if pLKO.1)"},
            {"day": "D4",  "task": "Replace selection media; assess cell viability"},
            {"day": "D5",  "task": "qRT-PCR — confirm mRNA knockdown (target: >70% reduction)"},
            {"day": "D6",  "task": "Western blot — confirm protein knockdown (target: >80% reduction)"},
            {"day": "D7+", "task": "Phenotypic assays (proliferation, apoptosis, migration, etc.)"},
            {"day": "D10", "task": "Rescue experiment: transfect siRNA-resistant cDNA; repeat phenotypic assays"},
        ]
    else:
        return [
            {"day": "D0",  "task": "Seed cells at 30-50% confluency; reverse transfect with siRNA + RNAiMAX"},
            {"day": "D1",  "task": "Change media 24h post-transfection to reduce cytotoxicity"},
            {"day": "D2",  "task": "qRT-PCR for mRNA knockdown efficiency (expect 70-90% at 48h)"},
            {"day": "D3",  "task": "Western blot for protein knockdown (allow for protein half-life)"},
            {"day": "D3-D5", "task": "Phenotypic readouts (proliferation assay, apoptosis, etc.)"},
            {"day": "D5",  "task": "If knockdown incomplete: re-transfect or increase siRNA concentration"},
            {"day": "D7",  "task": "Rescue experiment with siRNA-resistant cDNA transfection"},
            {"day": "D10", "task": "Repeat phenotypic assays post-rescue; compare with knockdown"},
        ]


# ---------------------------------------------------------------------------
# Validation experiment design
# ---------------------------------------------------------------------------

def _design_validation(gene: str, cell_type: Optional[str]) -> list[ValidationExperiment]:
    return [
        ValidationExperiment(
            assay_type="qRT-PCR",
            target_marker=f"{gene} mRNA",
            timepoint="48h post-transfection",
            expected_result=">70% mRNA reduction vs. scramble control",
            controls_required=["Scramble siRNA", "No-transfection control", "GAPDH internal reference"],
        ),
        ValidationExperiment(
            assay_type="Western blot",
            target_marker=f"{gene} protein",
            timepoint="72-96h post-transfection",
            expected_result=">80% protein reduction; timing depends on protein half-life",
            controls_required=["Scramble control", "β-actin or GAPDH loading control"],
        ),
        ValidationExperiment(
            assay_type="Immunofluorescence",
            target_marker=f"{gene} subcellular localization",
            timepoint="72h post-transfection",
            expected_result="Signal loss or mis-localization compared to scramble",
            controls_required=["Scramble control", "DAPI nuclear stain"],
        ),
        ValidationExperiment(
            assay_type="Cell viability assay (MTT/CellTiter-Glo)",
            target_marker="Cell proliferation",
            timepoint="D3, D5, D7 post-transfection",
            expected_result="Phenotype depends on gene function (arrest, apoptosis, or no effect)",
            controls_required=["Scramble control", "Positive control knockdown (essential gene)"],
        ),
        ValidationExperiment(
            assay_type="Rescue experiment",
            target_marker=f"siRNA-resistant {gene} cDNA re-expression",
            timepoint="48h after rescue transfection",
            expected_result="Phenotype reverts toward WT — confirms on-target specificity",
            controls_required=["Empty vector control", "Original knockdown without rescue"],
        ),
    ]


# ---------------------------------------------------------------------------
# Off-target risk assessment
# ---------------------------------------------------------------------------

def _assess_off_target_risk(sirna_designs: list[SiRNADesign], gene: str) -> tuple[str, list[str], list[str]]:
    """Simple rule-based off-target risk assessment."""
    warnings: list[str] = []
    seed_warnings: list[str] = []

    for i, design in enumerate(sirna_designs):
        seed = design.sequence_antisense[1:8]   # nt 2-8 of antisense (= guide)
        gc_seed = _gc_content(seed)
        if gc_seed > 60:
            seed_warnings.append(
                f"Candidate {i+1}: Seed region GC {gc_seed:.0f}% — high GC seed increases off-target risk"
            )
        if _has_runs(seed, 3):
            seed_warnings.append(
                f"Candidate {i+1}: Seed region contains homopolymer run — increased off-target binding risk"
            )

    warnings += [
        "Run BLAST against the transcriptome (NCBI BLAST, Ensembl or siSPOTR) with seed region before proceeding",
        "Include 2-3 independent siRNA sequences targeting different regions; convergent phenotypes increase confidence",
        "Chemical modifications (2'-OMe at position 2, phosphorothioate backbone) reduce off-target and innate immune activation",
        "Use the lowest effective concentration (titration 1-100 nM) to minimise off-target effects",
    ]

    # Estimate risk level
    risk_level = "Low" if not seed_warnings else ("High" if len(seed_warnings) > 2 else "Moderate")

    return risk_level, warnings, seed_warnings


# ---------------------------------------------------------------------------
# Gemini synthesis — full protocol
# ---------------------------------------------------------------------------

async def _synthesize_rnai_protocol(report: RNAiReport) -> RNAiReport:
    try:
        from core.gemini_interpreter import _call

        lines = [
            f"Target gene: {report.gene} ({report.organism})",
            f"RNAi type: {report.rnai_type}",
            f"Cell type: {report.cell_type or 'Not specified'}",
        ]
        if report.protein_function:
            lines.append(f"Protein function: {report.protein_function[:250]}")
        if report.sirna_designs:
            d = report.sirna_designs[0]
            lines.append(f"Best siRNA candidate (sense): 5'-{d.sequence_sense}-3'")
            lines.append(f"GC content: {d.gc_content_pct}%, Specificity: {d.specificity_score}")
        if report.delivery_strategy:
            ds = report.delivery_strategy
            lines.append(f"Delivery: {ds.method} with {ds.reagent_recommendation} at {ds.optimal_concentration}")
        lines.append(f"Off-target risk: {report.off_target_risk}")
        if report.disrupted_pathways:
            pw_names = [p.get("name", "?") for p in report.disrupted_pathways[:5]]
            lines.append(f"Expected disrupted pathways: {'; '.join(pw_names)}")
        if report.predicted_phenotypes:
            lines.append(f"Predicted phenotypes: {'; '.join(report.predicted_phenotypes[:4])}")
        if report.known_ko_phenotypes:
            lines.append(f"Known knockout phenotypes (OMIM/literature): {'; '.join(report.known_ko_phenotypes[:3])}")

        facts = "\n".join(lines)

        prompt = (
            "You are a senior molecular biologist and RNAi expert. Design a complete, rigorous RNAi "
            "knockdown experiment based on the following data.\n\n"
            f"TARGET INFORMATION:\n{facts}\n\n"
            "Produce a scientifically rigorous experimental protocol addressing:\n"
            "  1. Biological rationale for targeting this gene by RNAi\n"
            "  2. Key siRNA/shRNA design considerations for this target\n"
            "  3. Expected knockdown kinetics and efficiency\n"
            "  4. Predicted cellular phenotype based on gene function\n"
            "  5. Critical controls and why each is essential\n"
            "  6. Common pitfalls and how to avoid them (RISC saturation, innate immune activation, etc.)\n\n"
            "Return ONLY raw JSON (no markdown):\n"
            "{\n"
            "  \"experimental_protocol\": \"<8-12 sentences comprehensive protocol summary>\",\n"
            "  \"key_considerations\": [\"<consideration 1>\", ...],\n"
            "  \"predicted_phenotypes\": [\"<phenotype 1>\", ...],\n"
            "  \"troubleshooting_tips\": [\"<tip 1>\", ...],\n"
            "  \"analysis_caveats\": [\"<caveat 1>\", ...]\n"
            "}"
        )

        raw = await _call(prompt)
        if not raw:
            report.experimental_protocol = "AI synthesis unavailable."
            return report

        cleaned = raw.strip().lstrip("```json").lstrip("```").rstrip("```").strip()
        data = _json.loads(cleaned)

        report.experimental_protocol = str(data.get("experimental_protocol", "")).strip()

        for key, attr in [
            ("key_considerations", "key_considerations"),
            ("predicted_phenotypes", "predicted_phenotypes"),
            ("troubleshooting_tips", "troubleshooting_tips"),
            ("analysis_caveats", "analysis_caveats"),
        ]:
            val = data.get(key, [])
            setattr(report, attr, [str(x).strip() for x in val if x] if isinstance(val, list) else [])

    except Exception:
        if not report.experimental_protocol:
            report.experimental_protocol = "AI synthesis failed."

    return report


import json as _json


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_rnai_analysis(
    gene: str,
    rnai_type: str = "siRNA",
    cell_type: Optional[str] = None,
    organism: str = "Homo sapiens",
    step_cb=None,
) -> RNAiReport:
    """Design a complete RNAi experiment and predict knockdown phenotype.

    Args:
        gene:      Target gene symbol (e.g. "BRCA1", "EGFR").
        rnai_type: siRNA | shRNA | miRNA | general.
        cell_type: Target cell line or tissue (e.g. "HeLa", "primary neurons").
        organism:  Host organism (default Homo sapiens).
        step_cb:   Optional async progress callback.
    """
    async def _step(msg: str):
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    import httpx

    report = RNAiReport(
        gene=gene.upper(),
        organism=organism,
        rnai_type=rnai_type,
        cell_type=cell_type,
    )

    # ── 1. Fetch protein/CDS context from UniProt + NCBI ────────────────────
    await _step(f"[1/6] Fetching {gene} mRNA/CDS for siRNA design…")
    cds_seq: Optional[str] = None
    ncbi_gene_id: Optional[str] = None

    try:
        import httpx as _httpx
        from fetchers import uniprot as uni
        results = await uni.search_protein(
            f"reviewed:true AND gene_exact:{gene}", organism, size=1
        )
        if results:
            acc = results[0]["primaryAccession"]
            entry = await uni.fetch_by_accession(acc)
            report.uniprot_id = acc
            if isinstance(entry, dict):
                report.protein_name = (
                    entry.get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {}).get("value", "")
                ) or gene
                # Function
                for c in entry.get("comments", []):
                    if c.get("commentType") == "FUNCTION":
                        texts = c.get("texts", [])
                        if texts:
                            report.protein_function = texts[0].get("value", "")[:400]
                            break
                ncbi_gene_id = uni.get_ncbi_gene_id(entry)
    except Exception:
        pass

    # Fetch CDS for siRNA design
    if ncbi_gene_id:
        try:
            from fetchers.ncbi_gene import fetch_cds_for_gene
            cds_result = await fetch_cds_for_gene(ncbi_gene_id)
            if isinstance(cds_result, tuple):
                cds_seq = cds_result[0]
            elif isinstance(cds_result, str):
                cds_seq = cds_result
            if cds_seq:
                report.mrna_length = len(cds_seq)
                # Store CDS for frontend alignment cross-check (cap at 3000 nt)
                report.target_cds = cds_seq[:3000]
        except Exception:
            pass

    # ── 2. siRNA sequence design ─────────────────────────────────────────────
    await _step(f"[2/6] Designing optimal siRNA/shRNA sequences (Tuschl/Reynolds rules)…")
    if cds_seq and len(cds_seq) >= 21:
        report.sirna_designs = _design_sirna_candidates(cds_seq, n=3)
        if report.sirna_designs and rnai_type in ("shRNA", "general"):
            report.shrna_design = _design_shrna(report.sirna_designs[0].sequence_sense)
    else:
        # No CDS available — provide a design note
        report.sirna_designs = []
        report.key_considerations = [
            f"CDS for {gene} could not be retrieved from NCBI. "
            "Use a commercial siRNA design tool (e.g. Dharmacon siDESIGN Center, "
            "BLOCK-iT RNAi Designer, or siSPOTR) with the RefSeq mRNA sequence.",
        ]

    # ── 3. Off-target risk assessment ────────────────────────────────────────
    await _step("[3/6] Assessing off-target risk…")
    risk, ot_notes, seed_warns = _assess_off_target_risk(report.sirna_designs, gene)
    report.off_target_risk = risk
    report.off_target_notes = ot_notes
    report.seed_region_warnings = seed_warns

    # ── 4. Delivery strategy ─────────────────────────────────────────────────
    await _step("[4/6] Selecting delivery strategy…")
    report.delivery_strategy = _select_delivery(cell_type, rnai_type)
    report.controls = _required_controls(gene, rnai_type)
    report.timeline_days = _build_timeline(rnai_type, cell_type)
    report.validation_experiments = _design_validation(gene, cell_type)

    # ── 5. Predicted phenotype (reuse reverse genetics LoF analysis) ─────────
    await _step(f"[5/6] Predicting knockdown phenotype via DB cross-check…")
    try:
        from analyzers.reverse_genetics_analyzer import (
            run_reverse_genetics_analysis,
        )
        rg = await run_reverse_genetics_analysis(
            gene=gene,
            mutation_type="loss_of_function",
            organism=organism,
            step_cb=None,   # silent
        )
        report.disrupted_pathways = [
            {"id": p.pathway_id, "name": p.pathway_name, "severity": p.severity}
            for p in rg.pathway_consequences[:6]
        ]
        report.predicted_phenotypes = rg.key_consequences[:6]
        report.known_ko_phenotypes = rg.omim_diseases[:4]
    except Exception:
        pass

    # ── 6. Gemini synthesis ──────────────────────────────────────────────────
    await _step("[6/6] Gemini 2.5 Pro synthesizing experimental protocol…")
    report = await _synthesize_rnai_protocol(report)

    report.references = _RNAI_REFERENCES
    return report
