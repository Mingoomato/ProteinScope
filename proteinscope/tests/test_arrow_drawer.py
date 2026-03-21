"""Tests for annotators/arrow_drawer.py"""

from __future__ import annotations

import pytest
from PIL import Image, ImageDraw


def _make_draw(width: int = 400, height: int = 400):
    """Create a test PIL Image and Draw object."""
    img = Image.new("RGBA", (width, height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img)
    return img, draw


# ---------------------------------------------------------------------------
# draw_annotation
# ---------------------------------------------------------------------------

def test_draw_annotation_returns_bounding_box():
    """draw_annotation should return a 4-tuple bounding box."""
    from annotators.arrow_drawer import draw_annotation
    _, draw = _make_draw()
    box = draw_annotation(draw, 200, 200, 46, 17, 1, "Drug target: Warfarin", "top")
    assert isinstance(box, tuple)
    assert len(box) == 4
    x0, y0, x1, y1 = box
    assert x1 > x0
    assert y1 > y0


def test_draw_annotation_all_directions():
    """draw_annotation should not raise for any of the four directions."""
    from annotators.arrow_drawer import draw_annotation
    for direction in ["top", "bottom", "left", "right"]:
        _, draw = _make_draw()
        box = draw_annotation(draw, 200, 200, 46, 17, 1, "Test label", direction)
        assert len(box) == 4, f"Failed for direction: {direction}"


def test_draw_annotation_long_label_truncated():
    """Labels longer than MAX_LABEL_CHARS should be truncated with ellipsis."""
    from annotators.arrow_drawer import draw_annotation, MAX_LABEL_CHARS
    long_label = "A" * (MAX_LABEL_CHARS + 20)
    _, draw = _make_draw()
    # Should not raise; the truncation happens inside draw_annotation
    box = draw_annotation(draw, 200, 200, 46, 17, 1, long_label, "top")
    assert len(box) == 4


def test_draw_annotation_callout_numbers():
    """Multiple callout numbers should each be drawn without error."""
    from annotators.arrow_drawer import draw_annotation
    _, draw = _make_draw(800, 800)
    for i in range(1, 6):
        box = draw_annotation(draw, 100 + i * 100, 200, 46, 17, i, f"Label {i}", "top")
        assert len(box) == 4


# ---------------------------------------------------------------------------
# boxes_overlap
# ---------------------------------------------------------------------------

def test_boxes_overlap_true():
    """Two overlapping boxes should return True."""
    from annotators.arrow_drawer import boxes_overlap
    box_a = (0, 0, 100, 100)
    box_b = (50, 50, 150, 150)
    assert boxes_overlap(box_a, box_b) is True


def test_boxes_overlap_false():
    """Two non-overlapping boxes should return False."""
    from annotators.arrow_drawer import boxes_overlap
    box_a = (0, 0, 100, 100)
    box_b = (200, 200, 300, 300)
    assert boxes_overlap(box_a, box_b) is False


def test_boxes_overlap_touching_with_margin():
    """Boxes just outside margin should not overlap."""
    from annotators.arrow_drawer import boxes_overlap
    box_a = (0, 0, 100, 100)
    box_b = (110, 0, 210, 100)  # gap of 10 px > default margin 5
    assert boxes_overlap(box_a, box_b) is False


def test_boxes_overlap_within_margin():
    """Boxes within margin distance should be considered overlapping."""
    from annotators.arrow_drawer import boxes_overlap
    box_a = (0, 0, 100, 100)
    box_b = (103, 0, 203, 100)  # gap of 3 px < default margin 5
    assert boxes_overlap(box_a, box_b) is True


# ---------------------------------------------------------------------------
# _arrow_origin
# ---------------------------------------------------------------------------

def test_arrow_origin_top():
    """Arrow origin for 'top' should be above the node."""
    from annotators.arrow_drawer import _arrow_origin, ARROW_OFFSET
    ox, oy = _arrow_origin(100, 100, 46, 17, "top")
    assert ox == 100
    assert oy == 100 - 17 // 2 - ARROW_OFFSET


def test_arrow_origin_bottom():
    """Arrow origin for 'bottom' should be below the node."""
    from annotators.arrow_drawer import _arrow_origin, ARROW_OFFSET
    ox, oy = _arrow_origin(100, 100, 46, 17, "bottom")
    assert ox == 100
    assert oy == 100 + 17 // 2 + ARROW_OFFSET


def test_arrow_origin_left():
    """Arrow origin for 'left' should be to the left of the node."""
    from annotators.arrow_drawer import _arrow_origin, ARROW_OFFSET
    ox, oy = _arrow_origin(100, 100, 46, 17, "left")
    assert ox == 100 - 46 // 2 - ARROW_OFFSET
    assert oy == 100


def test_arrow_origin_right():
    """Arrow origin for 'right' should be to the right of the node."""
    from annotators.arrow_drawer import _arrow_origin, ARROW_OFFSET
    ox, oy = _arrow_origin(100, 100, 46, 17, "right")
    assert ox == 100 + 46 // 2 + ARROW_OFFSET
    assert oy == 100


def test_arrow_origin_unknown_defaults_to_top():
    """Unknown direction should fall back to 'top'."""
    from annotators.arrow_drawer import _arrow_origin, ARROW_OFFSET
    ox, oy = _arrow_origin(100, 100, 46, 17, "diagonal")
    assert ox == 100
    assert oy == 100 - 17 // 2 - ARROW_OFFSET


# ---------------------------------------------------------------------------
# _rounded_rect fallback
# ---------------------------------------------------------------------------

def test_rounded_rect_draws_without_error():
    """_rounded_rect should draw without raising even on old Pillow versions."""
    from annotators.arrow_drawer import _rounded_rect
    _, draw = _make_draw()
    # Should not raise regardless of Pillow version
    _rounded_rect(draw, [10, 10, 100, 50], radius=4,
                  fill=(255, 255, 200), outline=(200, 50, 50), width=1)


def test_rounded_rect_fallback_on_missing_method(monkeypatch):
    """_rounded_rect should fall back to rectangle if rounded_rectangle is missing."""
    from annotators import arrow_drawer
    from PIL import ImageDraw as _ImageDraw

    # Temporarily remove rounded_rectangle attribute if present
    original = getattr(_ImageDraw.ImageDraw, "rounded_rectangle", None)
    if original is not None:
        monkeypatch.delattr(_ImageDraw.ImageDraw, "rounded_rectangle", raising=False)

    _, draw = _make_draw()
    try:
        arrow_drawer._rounded_rect(draw, [10, 10, 100, 50], radius=4,
                                   fill=(255, 255, 200))
    except Exception as e:
        pytest.fail(f"_rounded_rect raised unexpectedly: {e}")
    finally:
        if original is not None:
            _ImageDraw.ImageDraw.rounded_rectangle = original


# ---------------------------------------------------------------------------
# _load_font
# ---------------------------------------------------------------------------

def test_load_font_returns_font_object():
    """_load_font should always return a font object (never raises)."""
    from annotators.arrow_drawer import _load_font
    font = _load_font("DejaVuSans.ttf", 11)
    # Just check it's a usable font object
    assert font is not None
