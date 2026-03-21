"""Master annotation pipeline — orchestrates the full annotation flow for one diagram.

Flow:
  1. Resolve node coordinates (KEGG KGML / Reactome JSON / PharmGKB SVG / STRING)
  2. Build rule-based comments per node (comment_engine)
  3. Polish comments with Gemini AI (ai_polisher) — concurrent
  4. Open image, draw arrows + callout bubbles (arrow_drawer)
  5. Save annotated PNG
  6. Build footnote block (footnote_builder)
  7. Return AnnotatedDiagram model
"""

from __future__ import annotations

import asyncio
import os
from pathlib import Path
from typing import Optional

from core.models import AnnotatedDiagram, NodeAnnotation, ProteinRecord

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

async def annotate_diagram(
    diagram_path: str,
    pathway_id: str,
    pathway_name: str,
    source: str,
    external_url: str,
    gene_symbols: list[str],
    record: ProteinRecord,
    save_dir: str,
    max_callouts: int = 10,
    no_ai_polish: bool = False,
) -> Optional[AnnotatedDiagram]:
    """Annotate a single pathway diagram PNG.

    Args:
        diagram_path:  Local path to the plain downloaded PNG.
        pathway_id:    Pathway identifier (KEGG map ID, Reactome stId, etc.).
        pathway_name:  Human-readable pathway name.
        source:        "KEGG" | "Reactome" | "PharmGKB" | "STRING".
        external_url:  Link to interactive viewer.
        gene_symbols:  Gene symbols to annotate on this diagram.
        record:        Full ProteinRecord for comment data.
        save_dir:      Directory to save the annotated PNG.
        max_callouts:  Maximum number of annotated nodes per diagram.
        no_ai_polish:  If True, skip AI polishing and use rule-based comments only.

    Returns:
        AnnotatedDiagram model, or None if annotation fails or no nodes found.
    """
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return None

    # ── Validate source image ────────────────────────────────────────────────
    if not diagram_path or not Path(diagram_path).exists():
        return None

    try:
        img = Image.open(diagram_path).convert("RGBA")
        img_w, img_h = img.size
    except Exception:
        return None

    if img_w < 50 or img_h < 50:
        return None

    # ── Step 1: Resolve node coordinates ────────────────────────────────────
    nodes = await _resolve_coordinates(
        source, pathway_id, gene_symbols, img_w, img_h
    )
    if not nodes:
        return None

    # Limit callout count
    nodes = nodes[:max_callouts]

    # ── Step 2: Build rule-based comments ───────────────────────────────────
    from annotators.callout_drawer import CalloutRegistry
    from annotators.comment_engine import (
        build_interaction_comment, build_drug_comment, build_disease_comment,
        build_stream_comment, build_arrow_label, get_protein_name,
    )

    registry = CalloutRegistry()
    raw_annotations: list[NodeAnnotation] = []

    for node in nodes:
        interaction = build_interaction_comment(node, record)
        drug        = build_drug_comment(node, record)
        disease     = build_disease_comment(node, record)
        up, down    = build_stream_comment(node, record)
        arrow_lbl   = build_arrow_label(interaction, drug)

        ann = NodeAnnotation(
            gene_symbol          = node["gene"],
            protein_name         = get_protein_name(node["gene"], record),
            x                    = node["x"],
            y                    = node["y"],
            node_width           = node["width"],
            node_height          = node["height"],
            arrow_label          = arrow_lbl,
            interaction_comment  = interaction,
            drug_comment         = drug,
            disease_comment      = disease,
            upstream_comment     = up,
            downstream_comment   = down,
            short_comment        = interaction,  # placeholder, replaced by AI below
            source               = source,
            diagram_source_url   = external_url,
        )
        registry.register(ann)
        raw_annotations.append(ann)

    # ── Step 3: AI polish (concurrent) ──────────────────────────────────────
    if not no_ai_polish:
        from annotators.ai_polisher import polish_comment
        tasks = [
            polish_comment({
                "gene_symbol":         ann.gene_symbol,
                "protein_name":        ann.protein_name,
                "interaction_comment": ann.interaction_comment,
                "drug_comment":        ann.drug_comment or "",
                "disease_comment":     ann.disease_comment or "",
                "upstream_comment":    ann.upstream_comment or "",
                "downstream_comment":  ann.downstream_comment or "",
            })
            for ann in raw_annotations
        ]
        polished_results = await asyncio.gather(*tasks, return_exceptions=True)
        for ann, result in zip(raw_annotations, polished_results):
            if isinstance(result, tuple) and len(result) == 2:
                ann.short_comment, ann.extended_footnote = result

    # ── Step 4: Draw annotations on image ───────────────────────────────────
    from annotators.arrow_drawer import draw_annotation, boxes_overlap

    draw = ImageDraw.Draw(img)
    occupied: list[tuple[int, int, int, int]] = []

    for ann in raw_annotations:
        direction = _choose_direction(ann.x, ann.y, img_w, img_h, occupied)
        try:
            box = draw_annotation(
                draw, ann.x, ann.y, ann.node_width, ann.node_height,
                ann.callout_number, ann.arrow_label, direction,
            )
            occupied.append(box)
        except Exception:
            pass  # Skip this annotation but continue with others

    # ── Step 5: Save annotated PNG ───────────────────────────────────────────
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    safe_id = pathway_id.replace(":", "_").replace("/", "_")
    out_filename = f"{source.lower()}_{safe_id}_annotated.png"
    out_path = str(Path(save_dir) / out_filename)

    try:
        img.convert("RGB").save(out_path, "PNG")
    except Exception:
        return None

    # ── Step 6: Build footnote ───────────────────────────────────────────────
    from annotators.footnote_builder import build_footnote_block
    footnote = build_footnote_block(raw_annotations, pathway_name, source)

    return AnnotatedDiagram(
        pathway_id           = pathway_id,
        pathway_name         = pathway_name,
        source               = source,
        plain_image_path     = diagram_path,
        annotated_image_path = out_path,
        external_url         = external_url,
        annotations          = raw_annotations,
        footnote_text        = footnote,
    )


# ---------------------------------------------------------------------------
# Direction chooser
# ---------------------------------------------------------------------------

def _choose_direction(
    x: int, y: int,
    img_w: int, img_h: int,
    occupied: list[tuple[int, int, int, int]],
) -> str:
    """Pick arrow direction that avoids existing label boxes and image edges."""
    from annotators.arrow_drawer import boxes_overlap, ARROW_OFFSET

    margin = 120
    # Priority order: top → bottom → right → left
    candidates = []
    if y > margin:             candidates.append("top")
    if y < img_h - margin:    candidates.append("bottom")
    if x < img_w - margin:    candidates.append("right")
    if x > margin:             candidates.append("left")
    if not candidates:
        candidates = ["top", "bottom", "right", "left"]

    for direction in candidates:
        # Approximate box for this direction
        if direction == "top":
            ox, oy = x, y - ARROW_OFFSET
        elif direction == "bottom":
            ox, oy = x, y + ARROW_OFFSET
        elif direction == "right":
            ox, oy = x + ARROW_OFFSET, y
        else:
            ox, oy = x - ARROW_OFFSET, y
        est_box = (ox - 80, oy - 15, ox + 80, oy + 15)
        if not any(boxes_overlap(est_box, occ) for occ in occupied):
            return direction

    return candidates[0]  # fallback: first candidate regardless of overlap


# ---------------------------------------------------------------------------
# Coordinate resolution dispatcher
# ---------------------------------------------------------------------------

async def _resolve_coordinates(
    source: str,
    pathway_id: str,
    gene_symbols: list[str],
    img_w: int,
    img_h: int,
) -> list[dict]:
    from annotators.coordinate_resolver import (
        get_kegg_node_coordinates,
        get_reactome_node_coordinates,
        get_pharmgkb_node_coordinates,
        get_string_node_coordinates,
    )
    try:
        if source == "KEGG":
            return await get_kegg_node_coordinates(
                pathway_id, gene_symbols, img_w, img_h)
        elif source == "Reactome":
            return await get_reactome_node_coordinates(
                pathway_id, gene_symbols, img_w, img_h)
        elif source == "PharmGKB":
            return await get_pharmgkb_node_coordinates(
                pathway_id, gene_symbols)
        elif source == "STRING":
            return await get_string_node_coordinates(
                gene_symbols[0] if gene_symbols else "",
                image_width=img_w, image_height=img_h)
        return []
    except Exception:
        return []
