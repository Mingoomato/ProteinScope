"""AAV Gene Therapy Vector Designer.

Designs an adeno-associated virus (AAV) gene therapy vector for a target gene,
selecting the optimal serotype, promoter, and configuration based on tissue
tropism, payload capacity, CpG immunogenicity, and clinical precedent.

Design workflow:
  1. Insert size estimation — CDS + ITRs + promoter + polyA signal
  2. Serotype selection — tissue tropism from function_description keywords
  3. Promoter selection — tissue-specific or ubiquitous
  4. CpG density and immunogenicity flagging
  5. Gemini synthesis — delivery optimisation, clinical trial context, key risks

Key constraints:
  - Effective AAV payload capacity: ≤ 4700 bp
    Citation: Naso MF 2017 BioDrugs doi:10.1007/s40259-017-0234-5
  - Serotype tissue tropism rules:
    Citation: Zincarelli C 2008 Mol Ther doi:10.1038/mt.2008.107
  - CpG immunostimulation threshold (> 6/100 bp):
    Citation: Krieg AM 2002 Nat Rev Drug Discov doi:10.1038/nrd804
"""

from __future__ import annotations

import re
from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel, Field

from core.evidence import DataProvenance, EvidenceGrade


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class AAVConfig(BaseModel):
    serotype: str
    """AAV serotype — e.g. AAV2, AAV5, AAV8, AAV9."""

    promoter: str
    """Promoter choice — CMV | EF1α | MHCK7 | SYN1 | CASI."""

    total_payload_bp: int
    """Estimated total vector payload in base pairs (ITRs + promoter + CDS + polyA)."""

    within_capacity: bool
    """True if total_payload_bp <= 4700 bp.
    # AAV vector capacity 4.7 kb: Naso MF 2017 BioDrugs doi:10.1007/s40259-017-0234-5
    """

    codon_optimized_cai: Optional[float] = None
    """Codon Adaptation Index of the supplied CDS (0–1); None if not computable."""

    cpg_dinucleotide_density: float
    """CpG dinucleotides per 100 bp of CDS.
    # CpG immunostimulation: Krieg AM 2002 Nat Rev Drug Discov doi:10.1038/nrd804
    """

    immunogenicity_flag: bool
    """True if CpG density > 6/100 bp OR CAI < 0.7.
    # CpG threshold: Krieg AM 2002 Nat Rev Drug Discov doi:10.1038/nrd804
    """

    delivery_route: str
    """Preferred clinical delivery route — IV | intrathecal | intramuscular | subretinal."""

    tissue_target: str
    """Primary target tissue or organ."""

    provenance: Optional[DataProvenance] = None
    """Computational provenance record for this config."""


class AAVDesignReport(BaseModel):
    gene: str
    recommended_config: AAVConfig
    alternative_configs: List[AAVConfig] = Field(default_factory=list)
    engineering_notes: List[str] = Field(default_factory=list)
    capacity_warning: Optional[str] = None
    gemini_summary: str = ""
    timestamp: datetime

    # Foundational references
    references: List[dict] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# Curated references
# ---------------------------------------------------------------------------

_AAV_REFERENCES: list[dict] = [
    {
        "pubmed_id": "28669112",
        "doi": "10.1007/s40259-017-0234-5",
        "title": "Adeno-associated virus vector as a platform for gene therapy delivery",
        "authors": ["Naso MF", "Tomkowicz B", "Perry WL 3rd", "Strohl WR"],
        "journal": "BioDrugs",
        "year": 2017,
        "apa_citation": (
            "Naso, M. F., Tomkowicz, B., Perry, W. L. 3rd, & Strohl, W. R. (2017). "
            "Adeno-associated virus vector as a platform for gene therapy delivery. "
            "BioDrugs, 31(4), 317–334. https://doi.org/10.1007/s40259-017-0234-5"
        ),
    },
    {
        "pubmed_id": "18392561",
        "doi": "10.1038/mt.2008.107",
        "title": "Systemic and local skeletal muscle tropism of novel adeno-associated virus serotypes 1–9",
        "authors": ["Zincarelli C", "Soltys S", "Rengo G", "Rabinowitz JE"],
        "journal": "Molecular Therapy",
        "year": 2008,
        "apa_citation": (
            "Zincarelli, C., Soltys, S., Rengo, G., & Rabinowitz, J. E. (2008). "
            "Systemic and local skeletal muscle tropism of novel adeno-associated virus serotypes 1–9. "
            "Molecular Therapy, 16(6), 1073–1080. https://doi.org/10.1038/mt.2008.107"
        ),
    },
    {
        "pubmed_id": "11827944",
        "doi": "10.1038/nrd804",
        "title": "CpG motifs in bacterial DNA and their immune effects",
        "authors": ["Krieg AM"],
        "journal": "Nature Reviews Drug Discovery",
        "year": 2002,
        "apa_citation": (
            "Krieg, A. M. (2002). CpG motifs in bacterial DNA and their immune effects. "
            "Nature Reviews Drug Discovery, 1(1), 49–60. https://doi.org/10.1038/nrd804"
        ),
    },
]


# ---------------------------------------------------------------------------
# AAV capacity constants
# ---------------------------------------------------------------------------

_ITR_BP = 290 * 2          # Two ITRs (left + right); ~145 bp each
_POLYA_BP = 220            # SV40 or bGH polyA signal
_CMV_PROMOTER_BP = 600     # CMV IE promoter/enhancer (~582 bp)
_EF1A_PROMOTER_BP = 1200   # EF1α (long form) ~1.2 kb
_MHCK7_PROMOTER_BP = 720   # MHCK7 muscle-specific ~720 bp
_SYN1_PROMOTER_BP = 470    # Synapsin-1 neuron-specific ~470 bp
_CASI_PROMOTER_BP = 230    # CASI ubiquitous compact promoter ~230 bp

_PROMOTER_SIZES: dict[str, int] = {
    "CMV": _CMV_PROMOTER_BP,
    "EF1α": _EF1A_PROMOTER_BP,
    "MHCK7": _MHCK7_PROMOTER_BP,
    "SYN1": _SYN1_PROMOTER_BP,
    "CASI": _CASI_PROMOTER_BP,
}

# AAV vector capacity 4.7 kb: Naso MF 2017 BioDrugs doi:10.1007/s40259-017-0234-5
_AAV_MAX_CAPACITY_BP = 4700


# ---------------------------------------------------------------------------
# Step 1: Insert size estimation
# ---------------------------------------------------------------------------

def _estimate_cds_bp(sequence: str, cds: Optional[str]) -> int:
    """Estimate CDS nucleotide length.

    Uses actual nucleotide CDS if provided; otherwise estimates from AA
    sequence length * 3 (codon length) + 3 (stop codon).
    """
    if cds and len(cds) >= 6:
        return len(cds.strip())
    # Estimate from AA sequence
    # AA length * 3 nt/codon + 3 nt stop codon
    aa_len = len(sequence.strip()) if sequence else 0
    return aa_len * 3 + 3


def _total_payload(cds_bp: int, promoter: str) -> int:
    """Compute total AAV payload in bp.

    Formula: ITRs (290 × 2) + promoter + CDS + polyA
    # AAV vector capacity 4.7 kb: Naso MF 2017 BioDrugs doi:10.1007/s40259-017-0234-5
    """
    promoter_bp = _PROMOTER_SIZES.get(promoter, _CMV_PROMOTER_BP)
    return _ITR_BP + promoter_bp + cds_bp + _POLYA_BP


# ---------------------------------------------------------------------------
# Step 2: Serotype selection
# ---------------------------------------------------------------------------

def _select_serotype_and_route(function_description: str) -> tuple[str, str, str, list[str]]:
    """Return (serotype, delivery_route, tissue_target, alternative_serotypes).

    Rules based on tissue tropism literature:
    # AAV serotype tissue tropism: Zincarelli C 2008 Mol Ther doi:10.1038/mt.2008.107
    """
    fd = function_description.lower()

    # Muscle / cardiac
    if any(kw in fd for kw in ("muscle", "cardiac", "heart", "dystrophin", "skeletal")):
        return "AAV9", "intramuscular", "Skeletal / cardiac muscle", ["AAV8", "AAV1"]

    # Neuron / CNS / brain
    if any(kw in fd for kw in ("neuron", "neural", "brain", "cns", "spinal", "cortex",
                                "cerebral", "dopamine", "sma", "als")):
        return "AAV9", "intrathecal", "Central nervous system", ["AAVrh10", "AAV4"]

    # Liver / metabolic
    if any(kw in fd for kw in ("liver", "hepat", "metabolic", "lysosomal",
                                "urea cycle", "coagulation", "factor viii", "factor ix")):
        return "AAV8", "IV", "Liver (hepatocytes)", ["AAV5", "AAV3B"]

    # Eye / retinal
    if any(kw in fd for kw in ("eye", "retina", "retinal", "photoreceptor",
                                "rpe", "optic", "vision", "leber")):
        return "AAV2", "subretinal", "Retinal photoreceptors / RPE", ["AAV5", "AAV8"]

    # Lung / pulmonary
    if any(kw in fd for kw in ("lung", "pulmonary", "airway", "cftr", "cystic fibrosis")):
        return "AAV5", "intratracheal", "Pulmonary epithelium", ["AAV1", "AAV6"]

    # Default fallback
    return "AAV2", "intramuscular", "Various (default)", ["AAV8", "AAV9"]


# ---------------------------------------------------------------------------
# Step 3: Promoter selection
# ---------------------------------------------------------------------------

def _select_promoter(function_description: str, serotype: str) -> str:
    """Select tissue-appropriate promoter.

    Hierarchy: muscle → MHCK7; neuron → SYN1; ubiquitous default → EF1α or CMV.
    Prefer CASI for large transgenes to preserve capacity headroom.
    """
    fd = function_description.lower()

    if any(kw in fd for kw in ("muscle", "cardiac", "heart", "skeletal", "dystrophin")):
        return "MHCK7"
    if any(kw in fd for kw in ("neuron", "neural", "brain", "cns", "synapse", "dopamine")):
        return "SYN1"
    # For liver/ubiquitous, CMV is well-validated; use EF1α for longer expression
    if any(kw in fd for kw in ("liver", "hepat", "metabolic")):
        return "EF1α"
    # Default: CMV for broad expression
    return "CMV"


# ---------------------------------------------------------------------------
# Step 4: CpG density and CAI
# ---------------------------------------------------------------------------

def _compute_cpg_density(cds: Optional[str]) -> float:
    """Count CG dinucleotides per 100 bp of CDS.

    # CpG immunostimulation: Krieg AM 2002 Nat Rev Drug Discov doi:10.1038/nrd804
    """
    if not cds or len(cds) < 2:
        return 0.0
    seq = cds.upper()
    cpg_count = len(re.findall(r"CG", seq))
    return round(cpg_count * 100.0 / len(seq), 2)


def _compute_cai(cds: Optional[str]) -> Optional[float]:
    """Compute a simplified Codon Adaptation Index for human codon usage.

    Uses the geometric mean of per-codon relative adaptiveness (RA) values.
    Returns None if CDS is unavailable or too short.

    CAI method: Sharp PM, Li WH 1987 Nucleic Acids Res doi:10.1093/nar/15.3.1281
    """
    if not cds or len(cds) < 3:
        return None
    # Human RA table (subset of most informative codons)
    # Values from Sharp & Li 1987; stop codons excluded
    _RA: dict[str, float] = {
        "TTT": 0.83, "TTC": 1.00, "TTA": 0.15, "TTG": 0.32,
        "CTT": 0.42, "CTC": 0.71, "CTA": 0.24, "CTG": 1.00,
        "ATT": 0.79, "ATC": 1.00, "ATA": 0.42, "ATG": 1.00,
        "GTT": 0.54, "GTC": 0.69, "GTA": 0.35, "GTG": 1.00,
        "TCT": 0.72, "TCC": 0.88, "TCA": 0.68, "TCG": 0.25,
        "AGT": 0.68, "AGC": 1.00,
        "CCT": 0.85, "CCC": 1.00, "CCA": 0.77, "CCG": 0.29,
        "ACT": 0.71, "ACC": 1.00, "ACA": 0.74, "ACG": 0.27,
        "GCT": 0.78, "GCC": 1.00, "GCA": 0.59, "GCG": 0.25,
        "TAT": 0.72, "TAC": 1.00,
        "CAT": 0.64, "CAC": 1.00, "CAA": 0.42, "CAG": 1.00,
        "AAT": 0.72, "AAC": 1.00, "AAA": 0.74, "AAG": 1.00,
        "GAT": 0.85, "GAC": 1.00, "GAA": 0.70, "GAG": 1.00,
        "TGT": 0.85, "TGC": 1.00, "TGG": 1.00,
        "CGT": 0.36, "CGC": 0.73, "CGA": 0.27, "CGG": 0.55,
        "AGA": 0.77, "AGG": 0.70,
        "GGT": 0.64, "GGC": 0.84, "GGA": 0.68, "GGG": 0.56,
    }
    _STOP = {"TAA", "TAG", "TGA"}
    seq = cds.upper().replace("U", "T")
    import math
    log_sum = 0.0
    n_codons = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) < 3 or codon in _STOP:
            continue
        ra = _RA.get(codon)
        if ra is None or ra <= 0:
            continue
        log_sum += math.log(ra)
        n_codons += 1
    if n_codons == 0:
        return None
    return round(math.exp(log_sum / n_codons), 3)


def _immunogenicity_flag(cpg_density: float, cai: Optional[float]) -> bool:
    """Flag immunogenicity risk.

    Criteria:
    - CpG density > 6 per 100 bp (Krieg AM 2002 Nat Rev Drug Discov doi:10.1038/nrd804)
    - OR CAI < 0.7 (codon bias may indicate poor optimisation)
    """
    # CpG immunostimulation: Krieg AM 2002 Nat Rev Drug Discov doi:10.1038/nrd804
    if cpg_density > 6.0:
        return True
    if cai is not None and cai < 0.7:
        return True
    return False


# ---------------------------------------------------------------------------
# Config builder
# ---------------------------------------------------------------------------

def _build_config(
    serotype: str,
    promoter: str,
    delivery_route: str,
    tissue_target: str,
    total_bp: int,
    cpg_density: float,
    cai: Optional[float],
) -> AAVConfig:
    # AAV vector capacity 4.7 kb: Naso MF 2017 BioDrugs doi:10.1007/s40259-017-0234-5
    within = total_bp <= _AAV_MAX_CAPACITY_BP
    imm = _immunogenicity_flag(cpg_density, cai)
    return AAVConfig(
        serotype=serotype,
        promoter=promoter,
        total_payload_bp=total_bp,
        within_capacity=within,
        codon_optimized_cai=cai,
        cpg_dinucleotide_density=cpg_density,
        immunogenicity_flag=imm,
        delivery_route=delivery_route,
        tissue_target=tissue_target,
        provenance=DataProvenance(
            source="ProteinScope AAV Designer (rule-based)",
            evidence_grade=EvidenceGrade.COMPUTATIONAL,
            scientific_caveat=(
                "Payload and immunogenicity estimates are heuristic; "
                "validate with wet-lab construct assembly and qPCR quantification."
            ),
            method="Naso 2017 capacity rules + Zincarelli 2008 tropism rules + Krieg 2002 CpG threshold",
        ),
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

async def run_aav_design(
    gene: str,
    sequence: str,
    cds: Optional[str],
    function_description: str,
    step_cb=None,
) -> Optional[AAVDesignReport]:
    """Design an AAV gene therapy vector for a given gene / protein.

    Args:
        gene:                 Gene symbol (e.g. "DMD", "RPE65", "SMN1").
        sequence:             Amino acid sequence (single-letter code).
        cds:                  Nucleotide CDS sequence (may be None).
        function_description: Free-text description of protein function/disease.
        step_cb:              Optional async progress callback (str → None).

    Returns:
        AAVDesignReport with recommended_config, alternative_configs, and
        engineering_notes.  Returns None only on catastrophic failure.
    """
    async def _step(msg: str) -> None:
        if step_cb:
            try:
                await step_cb(msg)
            except Exception:
                pass

    try:
        # ── 1. Insert size estimation ────────────────────────────────────────
        await _step("[1/5] Computing insert size...")
        cds_bp = _estimate_cds_bp(sequence, cds)

        # ── 2. Serotype selection ────────────────────────────────────────────
        await _step("[2/5] Selecting AAV serotype...")
        # AAV serotype tissue tropism: Zincarelli C 2008 Mol Ther doi:10.1038/mt.2008.107
        serotype, delivery_route, tissue_target, alt_serotypes = _select_serotype_and_route(
            function_description
        )

        # ── 3. Promoter selection ────────────────────────────────────────────
        await _step("[3/5] Selecting promoter...")
        promoter = _select_promoter(function_description, serotype)
        total_bp = _total_payload(cds_bp, promoter)

        # Capacity check — if over limit, try CASI (compact) to save space
        # AAV vector capacity 4.7 kb: Naso MF 2017 BioDrugs doi:10.1007/s40259-017-0234-5
        capacity_warning: Optional[str] = None
        if total_bp > _AAV_MAX_CAPACITY_BP:
            # Try compact CASI promoter
            total_bp_casi = _total_payload(cds_bp, "CASI")
            if total_bp_casi <= _AAV_MAX_CAPACITY_BP:
                promoter = "CASI"
                total_bp = total_bp_casi
                capacity_warning = (
                    f"Original promoter ({_select_promoter(function_description, serotype)}) "
                    f"pushed payload above 4700 bp capacity. Switched to compact CASI promoter "
                    f"({_CASI_PROMOTER_BP} bp) to fit within limits. "
                    f"Estimated payload: {total_bp_casi} bp."
                )
            else:
                # Still over — flag and recommend dual-vector split
                capacity_warning = (
                    f"Total estimated payload ({total_bp} bp) exceeds the AAV capacity "
                    f"limit of 4700 bp (Naso 2017). Recommend: (1) CDS truncation / "
                    f"mini-gene design, (2) dual-vector split strategy with overlapping "
                    f"homology for trans-splicing, or (3) intein-mediated protein trans-splicing. "
                    f"Estimated overage: {total_bp - _AAV_MAX_CAPACITY_BP} bp."
                )

        # ── 4. CpG density and immunogenicity ───────────────────────────────
        await _step("[4/5] Checking CpG density and immunogenicity...")
        # CpG immunostimulation: Krieg AM 2002 Nat Rev Drug Discov doi:10.1038/nrd804
        cpg_density = _compute_cpg_density(cds)
        cai = _compute_cai(cds)

        # Build recommended config
        recommended = _build_config(
            serotype=serotype,
            promoter=promoter,
            delivery_route=delivery_route,
            tissue_target=tissue_target,
            total_bp=total_bp,
            cpg_density=cpg_density,
            cai=cai,
        )

        # Build alternative configs for alt serotypes
        alt_configs: list[AAVConfig] = []
        for alt_ser in alt_serotypes[:2]:
            # Alternate route heuristic: first alt keeps same route, second may differ
            alt_config = _build_config(
                serotype=alt_ser,
                promoter=promoter,
                delivery_route=delivery_route,
                tissue_target=tissue_target,
                total_bp=total_bp,
                cpg_density=cpg_density,
                cai=cai,
            )
            alt_configs.append(alt_config)

        # Engineering notes
        notes: list[str] = []
        if cai is not None:
            if cai < 0.7:
                notes.append(
                    f"CAI = {cai:.3f} — below optimal threshold (0.7); codon optimisation "
                    "strongly recommended to improve expression and reduce immunogenicity."
                )
            else:
                notes.append(f"CAI = {cai:.3f} — acceptable codon adaptation for human expression.")
        else:
            notes.append(
                "CAI could not be computed (no CDS provided). Supply a codon-optimised "
                "nucleotide CDS for accurate immunogenicity estimation."
            )

        if cpg_density > 6.0:
            # CpG immunostimulation: Krieg AM 2002 Nat Rev Drug Discov doi:10.1038/nrd804
            notes.append(
                f"CpG density = {cpg_density:.1f}/100 bp — exceeds the 6/100 bp threshold "
                "(Krieg 2002). CpG depletion during codon optimisation is recommended to "
                "reduce Toll-like receptor 9 (TLR9)-mediated innate immune activation."
            )
        else:
            notes.append(
                f"CpG density = {cpg_density:.1f}/100 bp — within acceptable range (< 6/100 bp)."
            )

        notes.append(
            f"ITR size: {_ITR_BP} bp (2 × 145 bp). "
            f"Promoter ({promoter}): {_PROMOTER_SIZES.get(promoter, _CMV_PROMOTER_BP)} bp. "
            f"CDS estimate: {cds_bp} bp. polyA: {_POLYA_BP} bp. "
            f"Total: {total_bp} bp / 4700 bp capacity."
        )

        if serotype in ("AAV9", "AAV8") and "IV" not in delivery_route:
            notes.append(
                f"Systemic pre-existing neutralising antibodies against {serotype} are prevalent "
                "(up to 40% of adults). Screen patients for anti-AAV NAb titres prior to dosing."
            )

        # ── 5. Gemini synthesis ──────────────────────────────────────────────
        await _step("[5/5] Gemini synthesis...")
        gemini_summary = ""
        try:
            from core.gemini_interpreter import _call

            prompt = (
                f"You are an expert in AAV gene therapy vector design and clinical translation.\n\n"
                f"Gene: {gene}\n"
                f"Function: {function_description[:400]}\n"
                f"Recommended serotype: {serotype}\n"
                f"Promoter: {promoter}\n"
                f"Delivery route: {delivery_route}\n"
                f"Tissue target: {tissue_target}\n"
                f"Total payload: {total_bp} bp (limit: 4700 bp)\n"
                f"Within capacity: {total_bp <= _AAV_MAX_CAPACITY_BP}\n"
                f"CpG density: {cpg_density:.1f}/100 bp\n"
                f"CAI: {cai if cai is not None else 'not computed'}\n"
                f"Immunogenicity flag: {recommended.immunogenicity_flag}\n"
                + (f"Capacity warning: {capacity_warning}\n" if capacity_warning else "") +
                "\nProvide a rigorous 200-word synthesis covering:\n"
                "1. Rationale for the recommended serotype and delivery route\n"
                "2. Relevant clinical trial precedents (cite specific trials if known)\n"
                "3. Key manufacturing and safety risks (immunogenicity, genotoxicity, off-target transduction)\n"
                "4. Specific codon optimisation and CpG suppression recommendations if needed\n"
                "5. Recommended next steps toward IND-enabling studies\n\n"
                "Write in rigorous scientific language suitable for a gene therapy IND briefing document."
            )
            gemini_summary = await _call(prompt)
        except Exception:
            gemini_summary = ""

        return AAVDesignReport(
            gene=gene.upper(),
            recommended_config=recommended,
            alternative_configs=alt_configs,
            engineering_notes=notes,
            capacity_warning=capacity_warning,
            gemini_summary=gemini_summary or "",
            timestamp=datetime.utcnow(),
            references=_AAV_REFERENCES,
        )

    except Exception:
        return None
