"""Arrow drawer — Pillow-based arrow + callout bubble drawing on pathway PNG images.

All drawing functions degrade gracefully if fonts cannot be loaded.
"""

from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Optional

from PIL import Image, ImageDraw, ImageFont

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
ARROW_COLOR     = (220, 50, 50)       # red arrow and bubble
LABEL_BG_COLOR  = (255, 255, 200)     # pale yellow label background
LABEL_FG_COLOR  = (0, 0, 0)          # black label text
BUBBLE_COLOR    = (220, 50, 50)       # red filled circle
BUBBLE_TEXT_CLR = (255, 255, 255)     # white number in bubble

# ---------------------------------------------------------------------------
# Layout constants
# ---------------------------------------------------------------------------
ARROW_OFFSET    = 85    # px: how far from node edge the arrow origin sits
LABEL_PADDING   = 6     # px: padding inside label box
BUBBLE_RADIUS   = 13    # px: callout circle radius
MAX_LABEL_CHARS = 42    # truncate inline label if longer
FONT_SIZE_LABEL = 11
FONT_SIZE_BUBBLE = 11

# ---------------------------------------------------------------------------
# Font loading
# ---------------------------------------------------------------------------
_FONT_DIRS = [
    Path(__file__).parent.parent / "fonts",  # bundled: proteinscope/fonts/
    Path("fonts"),                            # relative fallback
]


def _load_font(filename: str, size: int) -> ImageFont.ImageFont:
    """Try to load a TTF font from known locations; fall back to PIL default."""
    for d in _FONT_DIRS:
        path = d / filename
        if path.exists():
            try:
                return ImageFont.truetype(str(path), size)
            except Exception:
                pass
    # System font fallback (Windows, macOS, Linux)
    system_paths = [
        "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/Arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for sp in system_paths:
        if os.path.exists(sp):
            try:
                return ImageFont.truetype(sp, size)
            except Exception:
                pass
    return ImageFont.load_default()


def _get_fonts() -> tuple[ImageFont.ImageFont, ImageFont.ImageFont]:
    label_font  = _load_font("DejaVuSans.ttf", FONT_SIZE_LABEL)
    bubble_font = _load_font("DejaVuSans-Bold.ttf", FONT_SIZE_BUBBLE)
    return label_font, bubble_font


# ---------------------------------------------------------------------------
# Public: draw single annotation
# ---------------------------------------------------------------------------

def draw_annotation(
    draw: ImageDraw.ImageDraw,
    node_x: int,
    node_y: int,
    node_w: int,
    node_h: int,
    callout_num: int,
    label_text: str,
    direction: str = "top",
) -> tuple[int, int, int, int]:
    """Draw numbered bubble + arrow + label box pointing at a node.

    Args:
        draw:         Pillow ImageDraw object.
        node_x/y:     Center coordinates of the target node.
        node_w/h:     Node dimensions (for computing edge offset).
        callout_num:  Integer shown in the bubble (1-based).
        label_text:   Inline short label text (truncated to MAX_LABEL_CHARS).
        direction:    Arrow direction: "top" | "bottom" | "left" | "right".

    Returns:
        (box_x0, box_y0, box_x1, box_y1) of the rendered label bounding box,
        for collision detection by the caller.
    """
    label_font, bubble_font = _get_fonts()

    # Compute arrow origin (outside node in given direction)
    origin_x, origin_y = _arrow_origin(node_x, node_y, node_w, node_h, direction)

    # Arrow line
    draw.line([(origin_x, origin_y), (node_x, node_y)],
              fill=ARROW_COLOR, width=2)

    # Arrowhead at node center
    _draw_arrowhead(draw, origin_x, origin_y, node_x, node_y, ARROW_COLOR)

    # Truncate label
    display_label = label_text[:MAX_LABEL_CHARS - 3] + "..." \
        if len(label_text) > MAX_LABEL_CHARS else label_text

    # Measure label text
    try:
        bbox = draw.textbbox((0, 0), display_label, font=label_font)
        text_w = bbox[2] - bbox[0]
        text_h = bbox[3] - bbox[1]
    except Exception:
        text_w, text_h = len(display_label) * 6, 14

    # Label box coordinates (centered at arrow origin)
    bx0 = origin_x - text_w // 2 - LABEL_PADDING
    by0 = origin_y - text_h // 2 - LABEL_PADDING
    bx1 = origin_x + text_w // 2 + LABEL_PADDING
    by1 = origin_y + text_h // 2 + LABEL_PADDING

    # Draw label background
    _rounded_rect(draw, [bx0, by0, bx1, by1], radius=4,
                  fill=LABEL_BG_COLOR, outline=ARROW_COLOR, width=1)

    # Draw label text
    try:
        draw.text((bx0 + LABEL_PADDING, by0 + LABEL_PADDING),
                  display_label, fill=LABEL_FG_COLOR, font=label_font)
    except Exception:
        draw.text((bx0 + LABEL_PADDING, by0 + LABEL_PADDING),
                  display_label, fill=LABEL_FG_COLOR)

    # Draw callout bubble with number (top-left of label box)
    bub_cx = bx0
    bub_cy = by0
    draw.ellipse(
        [bub_cx - BUBBLE_RADIUS, bub_cy - BUBBLE_RADIUS,
         bub_cx + BUBBLE_RADIUS, bub_cy + BUBBLE_RADIUS],
        fill=BUBBLE_COLOR,
    )
    num_str = str(callout_num)
    try:
        nbbox = draw.textbbox((0, 0), num_str, font=bubble_font)
        nw = nbbox[2] - nbbox[0]
        nh = nbbox[3] - nbbox[1]
        draw.text((bub_cx - nw // 2, bub_cy - nh // 2),
                  num_str, fill=BUBBLE_TEXT_CLR, font=bubble_font)
    except Exception:
        draw.text((bub_cx - 5, bub_cy - 7), num_str, fill=BUBBLE_TEXT_CLR)

    return bx0, by0, bx1, by1


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _arrow_origin(cx: int, cy: int, w: int, h: int,
                  direction: str) -> tuple[int, int]:
    """Compute arrow start point offset from node edge in the given direction."""
    half_w = w // 2
    half_h = h // 2
    offsets = {
        "top":    (cx,           cy - half_h - ARROW_OFFSET),
        "bottom": (cx,           cy + half_h + ARROW_OFFSET),
        "left":   (cx - half_w - ARROW_OFFSET, cy),
        "right":  (cx + half_w + ARROW_OFFSET, cy),
    }
    return offsets.get(direction, offsets["top"])


def _draw_arrowhead(draw: ImageDraw.ImageDraw,
                    x1: float, y1: float,
                    x2: float, y2: float,
                    color: tuple,
                    size: int = 8) -> None:
    """Draw a filled triangle arrowhead at (x2, y2) pointing from (x1, y1)."""
    angle = math.atan2(y2 - y1, x2 - x1)
    p1 = (x2 - size * math.cos(angle - math.pi / 6),
          y2 - size * math.sin(angle - math.pi / 6))
    p2 = (x2 - size * math.cos(angle + math.pi / 6),
          y2 - size * math.sin(angle + math.pi / 6))
    draw.polygon([(x2, y2), p1, p2], fill=color)


def _rounded_rect(draw: ImageDraw.ImageDraw, xy: list,
                  radius: int = 4, **kwargs) -> None:
    """Draw a rounded rectangle, falling back to regular rectangle if unavailable."""
    try:
        draw.rounded_rectangle(xy, radius=radius, **kwargs)
    except AttributeError:
        # Pillow < 8.2 does not have rounded_rectangle
        draw.rectangle(xy, **kwargs)


def boxes_overlap(box_a: tuple, box_b: tuple, margin: int = 5) -> bool:
    """Return True if two bounding boxes (x0,y0,x1,y1) overlap with margin."""
    ax0, ay0, ax1, ay1 = box_a
    bx0, by0, bx1, by1 = box_b
    return not (ax1 + margin < bx0 or bx1 + margin < ax0 or
                ay1 + margin < by0 or by1 + margin < ay0)
