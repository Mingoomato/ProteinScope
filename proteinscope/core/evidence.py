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
