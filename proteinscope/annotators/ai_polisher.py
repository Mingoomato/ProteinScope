"""AI polisher — use Gemini 2.5 Pro to produce polished annotation comments.

Uses the same Gemini pattern as core/gemini_interpreter.py.
Degrades gracefully: returns rule-based fallback if key is missing or API fails.
"""

from __future__ import annotations

import json as _json
from typing import Optional

THRESHOLD_CHARS = 300  # if total comment content exceeds this, produce extended footnote

SYSTEM_PROMPT = """You are a scientific annotator for a protein research report.
You will receive a JSON object describing a protein node in a biological pathway.
Your job is to write a concise, academically appropriate annotation.

RULES:
1. Write a SHORT COMMENT: 1-2 sentences max. Used as the footnote header.
2. Write an EXTENDED FOOTNOTE if total content warrants it (i.e. drug info,
   disease info, AND upstream/downstream are all present). The extended footnote
   is plain prose, 3-5 sentences, no bullet points.
3. Never invent data. Only use what is provided in the input JSON.
4. Use precise scientific language suitable for a thesis or research paper.
5. Return ONLY a JSON object with keys: "short_comment" and "extended_footnote".
   extended_footnote may be null if the short_comment covers everything.

FOOTNOTE TRIGGER THRESHOLD:
If the total character count of all comment fields combined exceeds 300 characters,
always produce an extended_footnote."""


async def polish_comment(node_data: dict) -> tuple[str, Optional[str]]:
    """Polish a rule-based node annotation using Gemini 2.5 Pro.

    Args:
        node_data: Dict with keys: gene_symbol, protein_name, interaction_comment,
                   drug_comment, disease_comment, upstream_comment, downstream_comment.

    Returns:
        (short_comment, extended_footnote | None)
        Falls back to rule-based interaction_comment on any failure.
    """
    total_chars = sum(
        len(str(v)) for v in node_data.values() if isinstance(v, str) and v
    )
    fallback_short = node_data.get("interaction_comment", "See pathway data.")

    try:
        # Lazy import to avoid circular dependency
        from core.gemini_interpreter import _call

        prompt = (
            f"{SYSTEM_PROMPT}\n\n"
            f"Input data:\n{_json.dumps(node_data, ensure_ascii=False)}\n\n"
            "Return ONLY the raw JSON object, no markdown fences:"
        )
        raw = await _call(prompt)
        if not raw:
            return _fallback(fallback_short, node_data, total_chars)

        # Strip markdown fences if present
        cleaned = raw.strip()
        if cleaned.startswith("```"):
            lines = cleaned.split("\n")
            cleaned = "\n".join(
                l for l in lines
                if not l.strip().startswith("```")
            ).strip()

        parsed = _json.loads(cleaned)
        short = str(parsed.get("short_comment") or fallback_short).strip()
        extended = parsed.get("extended_footnote")
        if extended:
            extended = str(extended).strip() or None
        return short, extended

    except Exception:
        return _fallback(fallback_short, node_data, total_chars)


def _fallback(
    short: str,
    node_data: dict,
    total_chars: int,
) -> tuple[str, Optional[str]]:
    """Return rule-based fallback when Gemini is unavailable."""
    if total_chars <= THRESHOLD_CHARS:
        return short, None
    parts = [
        node_data.get("drug_comment"),
        node_data.get("disease_comment"),
        node_data.get("upstream_comment"),
        node_data.get("downstream_comment"),
    ]
    extended = " ".join(p for p in parts if p) or None
    return short, extended
