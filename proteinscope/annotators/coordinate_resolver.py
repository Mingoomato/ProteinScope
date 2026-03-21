"""Coordinate resolver — extract pixel-level node positions from pathway diagram sources.

Supports: KEGG (KGML XML), Reactome (diagram layout JSON), STRING (network JSON API).
PharmGKB SVG parsing is included with a fallback to fixed position.
"""

from __future__ import annotations

import asyncio
from typing import Optional
import httpx
from lxml import etree


# ---------------------------------------------------------------------------
# KEGG KGML
# ---------------------------------------------------------------------------

async def get_kegg_node_coordinates(
    pathway_id: str,
    gene_symbols: list[str],
    image_width: int = 0,
    image_height: int = 0,
) -> list[dict]:
    """Extract node pixel coordinates from KEGG KGML XML.

    Args:
        pathway_id:   KEGG pathway ID, e.g. "hsa00010".
        gene_symbols: Gene symbols to locate in the KGML.
        image_width:  Actual downloaded PNG width (px). Used to compute scale.
        image_height: Actual downloaded PNG height (px).

    Returns:
        List of dicts with keys: gene, entry_id, x, y, width, height, source, pathway_id.
    """
    url = f"https://rest.kegg.jp/get/{pathway_id}/kgml"
    try:
        async with httpx.AsyncClient(timeout=30) as client:
            r = await client.get(url)
        if r.status_code != 200:
            return []
        root = etree.fromstring(r.content)
    except Exception:
        return []

    # Compute scale factor from KGML canvas size vs actual image size
    scale_x = scale_y = 1.0
    image_elem = root.find("image")
    if image_elem is not None and image_width > 0:
        kgml_w = float(image_elem.get("width", image_width) or image_width)
        kgml_h = float(image_elem.get("height", image_height) or image_height)
        if kgml_w > 0:
            scale_x = image_width / kgml_w
        if kgml_h > 0 and image_height > 0:
            scale_y = image_height / kgml_h

    results = []
    seen = set()

    for entry in root.findall("entry"):
        graphics = entry.find("graphics")
        if graphics is None:
            continue
        name = graphics.get("name", "")
        for gene in gene_symbols:
            if gene.upper() in name.upper() and gene not in seen:
                seen.add(gene)
                x = int(float(graphics.get("x", 0)) * scale_x)
                y = int(float(graphics.get("y", 0)) * scale_y)
                w = int(float(graphics.get("width", 46)) * scale_x)
                h = int(float(graphics.get("height", 17)) * scale_y)
                results.append({
                    "gene":       gene,
                    "entry_id":   entry.get("id", ""),
                    "x":          max(x, 10),
                    "y":          max(y, 10),
                    "width":      max(w, 20),
                    "height":     max(h, 10),
                    "source":     "KEGG",
                    "pathway_id": pathway_id,
                })
    return results


# ---------------------------------------------------------------------------
# Reactome diagram layout JSON
# ---------------------------------------------------------------------------

async def get_reactome_node_coordinates(
    pathway_id: str,
    gene_symbols: list[str],
    image_width: int = 0,
    image_height: int = 0,
) -> list[dict]:
    """Extract node pixel coordinates from Reactome diagram layout JSON.

    Args:
        pathway_id:   Reactome stable ID, e.g. "R-HSA-1640170".
        gene_symbols: Gene symbols to locate.
        image_width:  Actual downloaded PNG width (px).
        image_height: Actual downloaded PNG height (px).
    """
    url = f"https://reactome.org/ContentService/exporter/diagram/{pathway_id}.json"
    try:
        async with httpx.AsyncClient(timeout=30) as client:
            r = await client.get(url)
        if r.status_code != 200:
            return []
        layout = r.json()
    except Exception:
        return []

    # Compute scale from canvas bounds to image dimensions
    scale_x = scale_y = 1.0
    min_x = layout.get("minX", 0) or 0
    min_y = layout.get("minY", 0) or 0
    max_x = layout.get("maxX", 0) or 0
    max_y = layout.get("maxY", 0) or 0
    canvas_w = max_x - min_x
    canvas_h = max_y - min_y
    if canvas_w > 0 and image_width > 0:
        scale_x = image_width / canvas_w
    if canvas_h > 0 and image_height > 0:
        scale_y = image_height / canvas_h

    results = []
    seen = set()
    for node in layout.get("nodes", []):
        display_name = node.get("displayName", "")
        schema_class = node.get("schemaClass", "")
        if schema_class not in ("Protein", "EntityWithAccessionedSequence", "Complex", "Gene"):
            continue
        for gene in gene_symbols:
            if gene.upper() in display_name.upper() and gene not in seen:
                seen.add(gene)
                raw_x = float(node.get("x", 0) or 0) - min_x
                raw_y = float(node.get("y", 0) or 0) - min_y
                results.append({
                    "gene":       gene,
                    "entry_id":   node.get("stId", ""),
                    "x":          max(int(raw_x * scale_x), 10),
                    "y":          max(int(raw_y * scale_y), 10),
                    "width":      max(int(float(node.get("width", 60) or 60) * scale_x), 20),
                    "height":     max(int(float(node.get("height", 25) or 25) * scale_y), 10),
                    "source":     "Reactome",
                    "pathway_id": pathway_id,
                })
    return results


# ---------------------------------------------------------------------------
# PharmGKB SVG
# ---------------------------------------------------------------------------

async def get_pharmgkb_node_coordinates(
    pathway_id: str,
    gene_symbols: list[str],
) -> list[dict]:
    """Parse PharmGKB pathway SVG to extract node positions.

    Falls back to fixed placement if SVG parsing fails.
    """
    svg_url = f"https://api.pharmgkb.org/v1/data/pathway/{pathway_id}/download?format=svg"
    try:
        async with httpx.AsyncClient(timeout=60) as client:
            r = await client.get(svg_url)
        if r.status_code != 200:
            return _pharmgkb_fallback(pathway_id, gene_symbols)
        root = etree.fromstring(r.content)
    except Exception:
        return _pharmgkb_fallback(pathway_id, gene_symbols)

    results = []
    seen = set()

    # Walk all text elements looking for gene symbol matches
    for text_elem in root.iter("{http://www.w3.org/2000/svg}text"):
        content = (text_elem.text or "").strip()
        for gene in gene_symbols:
            if gene.upper() == content.upper() and gene not in seen:
                seen.add(gene)
                # Try to get position from element or parent group
                x = _svg_coord(text_elem, "x")
                y = _svg_coord(text_elem, "y")
                if x is None or y is None:
                    parent = text_elem.getparent()
                    if parent is not None:
                        transform = parent.get("transform", "")
                        x, y = _parse_svg_translate(transform)
                if x is not None and y is not None:
                    results.append({
                        "gene":       gene,
                        "entry_id":   "",
                        "x":          max(int(x), 10),
                        "y":          max(int(y), 10),
                        "width":      60,
                        "height":     20,
                        "source":     "PharmGKB",
                        "pathway_id": pathway_id,
                    })
    return results if results else _pharmgkb_fallback(pathway_id, gene_symbols)


def _svg_coord(elem, attr: str) -> Optional[float]:
    val = elem.get(attr)
    try:
        return float(val) if val else None
    except (ValueError, TypeError):
        return None


def _parse_svg_translate(transform: str) -> tuple[Optional[float], Optional[float]]:
    """Extract x,y from SVG transform='translate(x,y)' string."""
    import re
    m = re.search(r"translate\(\s*([\d.+-]+)[,\s]+([\d.+-]+)\s*\)", transform)
    if m:
        return float(m.group(1)), float(m.group(2))
    return None, None


def _pharmgkb_fallback(pathway_id: str, gene_symbols: list[str]) -> list[dict]:
    """Place annotations in a grid at top-right when SVG parsing fails."""
    results = []
    for i, gene in enumerate(gene_symbols[:5]):
        results.append({
            "gene":       gene,
            "entry_id":   "",
            "x":          400 + (i % 2) * 150,
            "y":          100 + (i // 2) * 80,
            "width":      60,
            "height":     20,
            "source":     "PharmGKB",
            "pathway_id": pathway_id,
        })
    return results


# ---------------------------------------------------------------------------
# STRING network JSON
# ---------------------------------------------------------------------------

async def get_string_node_coordinates(
    gene_symbol: str,
    species: int = 9606,
    image_width: int = 1000,
    image_height: int = 1000,
) -> list[dict]:
    """Get STRING network node positions via the STRING API.

    The STRING network PNG default size is 1000x1000 px, and the API
    returns normalized coordinates (0-1 range or absolute, depending on version).
    We request the network layout and scale to image dimensions.
    """
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": gene_symbol,
        "species": species,
        "network_type": "functional",
        "required_score": 400,
    }
    try:
        async with httpx.AsyncClient(timeout=30) as client:
            r = await client.get(url, params=params)
        if r.status_code != 200:
            return []
        data = r.json()
        if not isinstance(data, list):
            return []
    except Exception:
        return []

    results = []
    seen = set()
    for node in data:
        name = node.get("preferredName", node.get("stringId", ""))
        if name.upper() == gene_symbol.upper() and gene_symbol not in seen:
            seen.add(gene_symbol)
            # STRING API returns x,y in the range corresponding to network layout
            x = float(node.get("x", image_width / 2) or image_width / 2)
            y = float(node.get("y", image_height / 2) or image_height / 2)
            results.append({
                "gene":       gene_symbol,
                "entry_id":   node.get("stringId", ""),
                "x":          max(int(x), 10),
                "y":          max(int(y), 10),
                "width":      40,
                "height":     40,
                "source":     "STRING",
                "pathway_id": f"string_{gene_symbol}",
            })
    return results
