"""ESM-2 masked marginal likelihood variant scorer.

Scores missense variants using Facebook's ESM-2 protein language model.
The masked marginal likelihood score (ΔLL) measures how much each
substitution reduces the model's confidence in the wild-type residue —
a proxy for functional impact.

CUDA acceleration
-----------------
When a CUDA-capable GPU is detected via `utils.device.get_device()`, the
ESM-2 model is loaded onto GPU automatically.  On CPU-only hosts the same
code runs correctly, just slower.

Installation
------------
ESM-2 requires the ``fair-esm`` package::

    pip install fair-esm

If fair-esm is not installed, all functions return gracefully with None /
empty results.  No ImportError is raised.

Usage::

    from analyzers.variant_scorer import score_variants

    results = score_variants(
        sequence="MKTAY...",
        variants=[("A", 42, "V"), ("G", 100, "D")],
    )
    for r in results:
        print(r.variant, r.delta_log_likelihood, r.predicted_impact)
"""

from __future__ import annotations

import logging
from typing import Optional

from pydantic import BaseModel, Field

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Result model
# ---------------------------------------------------------------------------


class VariantScore(BaseModel):
    """Scored missense variant."""

    position: int          # 1-based residue position
    ref_aa: str            # wild-type amino acid (single-letter)
    alt_aa: str            # substituted amino acid (single-letter)
    variant: str           # e.g. "A42V"
    delta_log_likelihood: Optional[float] = None  # ΔLL; more negative = more damaging
    wild_type_ll: Optional[float] = None  # LL of WT residue at this position
    alt_ll: Optional[float] = None        # LL of mutant residue at this position
    predicted_impact: Optional[str] = None  # "Benign" / "Possibly damaging" / "Likely damaging"
    model_name: str = ""   # ESM model used (e.g. "esm2_t33_650M_UR50D")


# ---------------------------------------------------------------------------
# Internal: model cache
# ---------------------------------------------------------------------------

_esm_cache: dict = {}  # {model_name: (model, alphabet, batch_converter)}


def _get_esm_model(model_name: str = "esm2_t33_650M_UR50D"):
    """Load and cache an ESM-2 model.  Returns (model, alphabet, batch_converter) or None."""
    if model_name in _esm_cache:
        return _esm_cache[model_name]

    try:
        import esm
        from utils.device import get_device

        device = get_device()
        _log.info("Loading ESM-2 model '%s' on device=%s …", model_name, device)

        model, alphabet = esm.pretrained.__dict__[model_name]()
        model = model.eval().to(device)
        batch_converter = alphabet.get_batch_converter()

        _esm_cache[model_name] = (model, alphabet, batch_converter)
        _log.info("ESM-2 '%s' loaded successfully", model_name)
        return _esm_cache[model_name]
    except (ImportError, KeyError, Exception) as exc:
        _log.warning("ESM-2 model load failed (%s): %s", model_name, exc)
        return None


# ---------------------------------------------------------------------------
# Scoring logic
# ---------------------------------------------------------------------------


def _impact_label(delta_ll: float) -> str:
    """Map ΔLL to an interpretable impact category.

    Thresholds calibrated on ClinVar/ProteinGym benchmarks (approximate):
        ΔLL > -1   → Benign
        -4 < ΔLL ≤ -1 → Possibly damaging
        ΔLL ≤ -4   → Likely damaging
    """
    if delta_ll > -1.0:
        return "Benign"
    if delta_ll > -4.0:
        return "Possibly damaging"
    return "Likely damaging"


def score_variants(
    sequence: str,
    variants: list[tuple[str, int, str]],
    model_name: str = "esm2_t33_650M_UR50D",
) -> list[VariantScore]:
    """Score a list of missense variants using ESM-2 masked marginal likelihood.

    Args:
        sequence:   Wild-type protein sequence (single-letter amino acids).
        variants:   List of (ref_aa, position_1based, alt_aa) tuples.
                    e.g. [("A", 42, "V"), ("G", 100, "D")]
        model_name: ESM-2 model to use.  Smaller models are faster:
                    - esm2_t6_8M_UR50D    (8 M params — fastest)
                    - esm2_t12_35M_UR50D  (35 M params)
                    - esm2_t33_650M_UR50D (650 M params — default, best accuracy)
                    - esm2_t36_3B_UR50D   (3 B params — requires high VRAM)

    Returns:
        List of VariantScore objects, one per variant.
        Returns empty list if fair-esm is not installed or model load fails.
    """
    if not variants or not sequence:
        return []

    esm_pkg = _get_esm_model(model_name)
    if esm_pkg is None:
        return [
            VariantScore(
                position=pos,
                ref_aa=ref,
                alt_aa=alt,
                variant=f"{ref}{pos}{alt}",
                model_name=model_name,
            )
            for ref, pos, alt in variants
        ]

    try:
        import torch
        from utils.device import get_device

        device = get_device()
        model, alphabet, batch_converter = esm_pkg

        results: list[VariantScore] = []

        with torch.no_grad():
            for ref_aa, position, alt_aa in variants:
                pos_0 = position - 1  # convert to 0-based
                if pos_0 < 0 or pos_0 >= len(sequence):
                    results.append(VariantScore(
                        position=position, ref_aa=ref_aa, alt_aa=alt_aa,
                        variant=f"{ref_aa}{position}{alt_aa}", model_name=model_name,
                    ))
                    continue

                # Mask the target position
                masked_seq = list(sequence)
                masked_seq[pos_0] = "<mask>"
                masked_str = "".join(masked_seq)

                _, _, batch_tokens = batch_converter([("protein", masked_str)])
                batch_tokens = batch_tokens.to(device)

                output = model(batch_tokens, repr_layers=[], return_contacts=False)
                logits = output["logits"]  # shape: [1, L+2, vocab_size]

                # +1 offset because ESM prepends a <cls> token
                tok_idx = pos_0 + 1
                log_probs = torch.log_softmax(logits[0, tok_idx], dim=-1)

                ref_idx = alphabet.get_idx(ref_aa)
                alt_idx = alphabet.get_idx(alt_aa)

                wt_ll = float(log_probs[ref_idx].item())
                alt_ll = float(log_probs[alt_idx].item())
                delta_ll = alt_ll - wt_ll

                results.append(VariantScore(
                    position=position,
                    ref_aa=ref_aa,
                    alt_aa=alt_aa,
                    variant=f"{ref_aa}{position}{alt_aa}",
                    delta_log_likelihood=round(delta_ll, 4),
                    wild_type_ll=round(wt_ll, 4),
                    alt_ll=round(alt_ll, 4),
                    predicted_impact=_impact_label(delta_ll),
                    model_name=model_name,
                ))

        return results

    except Exception as exc:
        _log.warning("ESM-2 variant scoring failed: %s", exc)
        return []


def score_all_possible_variants(
    sequence: str,
    positions: list[int] | None = None,
    model_name: str = "esm2_t33_650M_UR50D",
) -> list[VariantScore]:
    """Score every possible single-amino-acid substitution at given positions.

    If positions is None, scores all positions in the sequence (very slow for
    long sequences — use on GPU or limit to active site residues).

    Args:
        sequence:   Wild-type protein sequence.
        positions:  1-based positions to scan.  None → all positions.
        model_name: ESM-2 model.

    Returns:
        List of VariantScore for all (position, alt_aa) combinations.
    """
    AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

    if positions is None:
        positions = list(range(1, len(sequence) + 1))

    variants: list[tuple[str, int, str]] = []
    for pos in positions:
        ref = sequence[pos - 1]
        for alt in AA_ALPHABET:
            if alt != ref:
                variants.append((ref, pos, alt))

    return score_variants(sequence, variants, model_name=model_name)
