"""Evidence grading and data provenance for the Feature Confidence Framework.

Every data point surfaced by ProteinScope carries an EvidenceGrade and
optionally a DataProvenance record so the user always knows *how confident*
a result is and *where it came from*.

Import pattern (in analyzers and models)::

    from core.evidence import EvidenceGrade, DataProvenance
"""

from __future__ import annotations

from datetime import datetime
from enum import Enum
from typing import Optional

from pydantic import BaseModel


class EvidenceGrade(str, Enum):
    """Confidence tier for a single data point.

    Listed from highest to lowest evidential strength.
    """

    EXPERIMENTAL = "experimental"
    """Derived from direct wet-lab measurement.

    Examples: X-ray / cryo-EM crystal structure, in-vitro binding assay,
    ClinVar pathogenicity assertion with functional evidence, published
    biochemical characterisation.
    """

    COMPUTATIONAL = "computational"
    """Derived from a computational prediction or model.

    Examples: AlphaFold pLDDT / 3-D coordinates, ESM-2 ΔLL variant score,
    PROPKA pKa prediction, IUPred3 disorder score, FuzDrop LLPS propensity.
    """

    LITERATURE = "literature"
    """Extracted or inferred from published literature by an LLM.

    Examples: gene-disease associations parsed from PubMed abstracts,
    Gemini-extracted mechanism of action, pathway membership derived from
    review articles.
    """

    AI_GENERATED = "ai_generated"
    """Synthesised de-novo by a generative model (Gemini) with no direct
    literature anchor.

    Examples: Gemini 'therapeutic_implications' narrative, AI-generated
    experimental recommendation, Gemini regulatory summary.
    """


class DataProvenance(BaseModel):
    """Provenance record attached to a single data point or analysis result.

    Attach one of these to any model field where the user should understand
    the origin, confidence, and caveats of the data.

    Example::

        provenance = DataProvenance(
            source="UniProt REST v2",
            evidence_grade=EvidenceGrade.EXPERIMENTAL,
            scientific_caveat=None,
        )

        provenance = DataProvenance(
            source="AlphaFold EBI v4",
            evidence_grade=EvidenceGrade.COMPUTATIONAL,
            scientific_caveat="AlphaFold prediction; no crystal structure deposited",
        )

        provenance = DataProvenance(
            source="ESM-2 650M (fair-esm)",
            evidence_grade=EvidenceGrade.COMPUTATIONAL,
            confidence_interval="ΔLL calibration ±0.3 on ProteinGym benchmark",
            scientific_caveat="Zero-shot prediction; not validated for this specific variant",
        )
    """

    source: str
    """Human-readable name of the data source.

    Examples: "UniProt REST v2", "ESM-2 650M", "AlphaFold EBI v4",
    "PhosphoSitePlus GOLD", "IUPred3", "Gemini 2.5 Pro synthesis".
    """

    retrieved_at: Optional[datetime] = None
    """UTC timestamp when the data was fetched (populated by fetchers)."""

    evidence_grade: EvidenceGrade
    """Confidence tier — see EvidenceGrade for definitions."""

    confidence_interval: Optional[str] = None
    """Free-text CI / error estimate.

    Examples: "±2 kcal/mol (MM/GBSA)", "pLDDT 91.3 (high confidence)",
    "ΔLL calibration ±0.3 (ProteinGym)".
    """

    scientific_caveat: Optional[str] = None
    """One-sentence caveat the user should know before acting on this data.

    Examples:
    - "AlphaFold prediction; no crystal structure deposited in PDB"
    - "ESM-2 zero-shot — not validated for this specific variant"
    - "IUPred3 per-residue score; experimental IDR boundaries may differ"
    """

    method: Optional[str] = None
    """Specific computational method or assay type, if applicable.

    Examples: "BLOSUM62 global alignment", "masked marginal likelihood",
    "Elastic Network Model (ANM)", "Western blot".
    """


class QBConfidenceScore(BaseModel):
    """Quantum Biology Confidence Score — 윤박사 (CALTECH), Grand Consortium V3 2026-03-20.

    Formula: score = (min_pLDDT_path / 100) * cofactor_certainty * kie_evidence_flag

    Uses MIN pLDDT across all ETP pathway residues (not average). A single
    low-confidence residue makes the entire exponential decay calculation
    unreliable because k_ET ∝ exp(-βR) is hypersensitive to distance errors.

    cofactor_certainty tiers:
      1.0 — experimentally confirmed in this protein (PDB / biochemical assay)
      0.6 — confirmed in close homolog (>90% identity)
      0.2 — predicted from sequence (domain annotation)
      0.0 — no evidence

    kie_evidence_flag tiers (kH/kD threshold kH/kD > 10 at 25°C):
      1.0 — kH/kD > 10 published for this enzyme
      0.5 — kH/kD 7–10, or > 10 in close homolog
      0.1 — no KIE data available

    Badge colours (윤박사):
      0.0–0.3 → red   "Speculative"       (insufficient structural/experimental support)
      0.3–0.7 → yellow "Moderate Evidence" (some data, additional validation needed)
      0.7–1.0 → green "Strong Evidence"   (reliable structure + cofactor + KIE)

    # Citation: Rodgers & Hore (2009) PNAS 106:353-360. doi:10.1073/pnas.0711968106
    # Citation: Cha, Murray & Klinman (1989) Science 243:1325. doi:10.1126/science.243.4896.1325
    """

    min_plddt_path: float          # minimum pLDDT across all ETP pathway residues
    cofactor_certainty: float      # 1.0 | 0.6 | 0.2 | 0.0
    kie_evidence_flag: float       # 1.0 | 0.5 | 0.1
    score: float                   # = (min_plddt_path/100) * cofactor_certainty * kie_evidence_flag
    badge_color: str               # "red" | "yellow" | "green"
    badge_label: str               # "Speculative" | "Moderate Evidence" | "Strong Evidence"

    @classmethod
    def compute(
        cls,
        min_plddt: float,
        cofactor_certainty: float,
        kie_flag: float,
    ) -> "QBConfidenceScore":
        """Compute score and assign badge from raw inputs."""
        score = (min_plddt / 100.0) * cofactor_certainty * kie_flag
        score = round(min(max(score, 0.0), 1.0), 3)
        if score >= 0.7:
            color, label = "green", "Strong Evidence"
        elif score >= 0.3:
            color, label = "yellow", "Moderate Evidence"
        else:
            color, label = "red", "Speculative"
        return cls(
            min_plddt_path=min_plddt,
            cofactor_certainty=cofactor_certainty,
            kie_evidence_flag=kie_flag,
            score=score,
            badge_color=color,
            badge_label=label,
        )
