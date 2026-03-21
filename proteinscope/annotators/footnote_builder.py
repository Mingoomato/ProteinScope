"""Footnote builder — assemble the annotation legend below each pathway diagram."""

from __future__ import annotations

from core.models import NodeAnnotation

# Unicode circled numbers ①–⑮ for callout display
_CIRCLE_NUMS = "①②③④⑤⑥⑦⑧⑨⑩⑪⑫⑬⑭⑮"


def build_footnote_block(
    annotations: list[NodeAnnotation],
    pathway_name: str,
    source: str,
) -> str:
    """Build a formatted footnote legend block for Markdown reports.

    Format per annotation:
        ① GENE — Protein Name
           Short comment.
           ┗ Extended footnote sentence 1.
           ┗ Extended footnote sentence 2.

    Args:
        annotations:  List of NodeAnnotation objects in callout number order.
        pathway_name: Display name of the pathway.
        source:       Data source string (e.g. "KEGG", "Reactome").

    Returns:
        Formatted markdown string for the footnote block.
    """
    if not annotations:
        return ""

    lines: list[str] = [
        f"\n**Annotation Notes — {pathway_name} ({source})**\n"
    ]

    for ann in annotations:
        n = ann.callout_number
        num_char = _CIRCLE_NUMS[n - 1] if 1 <= n <= len(_CIRCLE_NUMS) else f"({n})"

        lines.append(f"{num_char} **{ann.gene_symbol}** — {ann.protein_name}")

        short = ann.short_comment or ann.interaction_comment
        if short:
            lines.append(f"   {short}")

        if ann.extended_footnote:
            # Split extended footnote into sentences and indent each
            sentences = [s.strip() for s in ann.extended_footnote.split(". ") if s.strip()]
            for sentence in sentences:
                if not sentence.endswith("."):
                    sentence += "."
                lines.append(f"   ┗ {sentence}")

        lines.append("")  # blank line between entries

    return "\n".join(lines)


def build_pdf_footnote_rows(annotations: list[NodeAnnotation]) -> list[dict]:
    """Build a list of dicts for PDF rendering.

    Each dict has: callout_num, gene_symbol, protein_name, short_comment, extended_footnote.
    """
    rows = []
    for ann in annotations:
        rows.append({
            "callout_num":       ann.callout_number,
            "gene_symbol":       ann.gene_symbol,
            "protein_name":      ann.protein_name,
            "short_comment":     ann.short_comment or ann.interaction_comment,
            "extended_footnote": ann.extended_footnote,
        })
    return rows
