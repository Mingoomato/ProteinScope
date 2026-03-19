"""Diagram fetcher — concurrent coordinator for all pathway image downloads.

Saves images locally and optionally annotates them with the annotation pipeline.
Returns both local path (for PDF embedding) and external URL (for clickable links).
"""

from __future__ import annotations

import asyncio
import os
from pathlib import Path
from typing import Optional


async def fetch_all_diagrams(
    record,
    save_dir: str = "./diagrams",
    annotate: bool = True,
    max_callouts: int = 10,
    no_ai_polish: bool = False,
) -> dict:
    """Download all pathway diagram images for a ProteinRecord.

    If annotate=True, each downloaded PNG is passed through the annotation
    pipeline and the record's annotated_diagrams list is populated.

    Args:
        record:        Assembled ProteinRecord.
        save_dir:      Root directory for images.
        annotate:      Whether to run the annotation pipeline after download.
        max_callouts:  Max annotated nodes per diagram.
        no_ai_polish:  Skip Gemini AI polishing (use rule-based comments only).

    Returns:
        Dict mapping source name to list of result dicts.
    """
    protein_dir = str(Path(save_dir) / record.uniprot_id)
    Path(protein_dir).mkdir(parents=True, exist_ok=True)

    annotate_dir = str(Path(protein_dir) / "annotated")

    tasks = {
        "pharmgkb": _download_pharmgkb_images(
            record, protein_dir,
            annotate=annotate, max_callouts=max_callouts,
            no_ai_polish=no_ai_polish, annotate_dir=annotate_dir,
        ),
        "kegg": _download_kegg_images(
            record, protein_dir,
            annotate=annotate, max_callouts=max_callouts,
            no_ai_polish=no_ai_polish, annotate_dir=annotate_dir,
        ),
        "reactome": _download_reactome_images(
            record, protein_dir,
            annotate=annotate, max_callouts=max_callouts,
            no_ai_polish=no_ai_polish, annotate_dir=annotate_dir,
        ),
        "string": _download_string_image(
            record, protein_dir,
            annotate=annotate, max_callouts=max_callouts,
            no_ai_polish=no_ai_polish, annotate_dir=annotate_dir,
        ),
    }

    results_list = await asyncio.gather(*tasks.values(), return_exceptions=True)
    results = {}
    for key, result in zip(tasks.keys(), results_list):
        results[key] = result if not isinstance(result, Exception) else []
    return results


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

async def _maybe_annotate(
    local_path: str,
    pathway_id: str,
    pathway_name: str,
    source: str,
    external_url: str,
    gene_symbols: list[str],
    record,
    annotate: bool,
    max_callouts: int,
    no_ai_polish: bool,
    annotate_dir: str,
) -> Optional[str]:
    """Run the annotation pipeline on a downloaded image.

    Updates record.annotated_diagrams and returns the annotated image path.
    Returns the original path if annotation is disabled or fails.
    """
    if not annotate:
        return local_path
    try:
        from annotators.annotation_pipeline import annotate_diagram
        result = await annotate_diagram(
            diagram_path = local_path,
            pathway_id   = pathway_id,
            pathway_name = pathway_name,
            source       = source,
            external_url = external_url,
            gene_symbols = gene_symbols,
            record       = record,
            save_dir     = annotate_dir,
            max_callouts = max_callouts,
            no_ai_polish = no_ai_polish,
        )
        if result:
            record.annotated_diagrams.append(result)
            return result.annotated_image_path
    except Exception:
        pass
    return local_path


async def _download_pharmgkb_images(
    record,
    save_dir: str,
    annotate: bool,
    max_callouts: int,
    no_ai_polish: bool,
    annotate_dir: str,
) -> list[dict]:
    from fetchers.pharmgkb import download_pathway_image
    gene_symbols = [record.gene_name]
    items = []
    for pathway in getattr(record, "pharmacogenomic_pathways", []):
        local_path = await download_pathway_image(pathway.pathway_id, save_dir)
        if local_path:
            final_path = await _maybe_annotate(
                local_path, pathway.pathway_id, pathway.pathway_name,
                "PharmGKB", pathway.diagram_url, gene_symbols,
                record, annotate, max_callouts, no_ai_polish, annotate_dir,
            )
            pathway.diagram_image_path = final_path
            items.append({
                "pathway_id": pathway.pathway_id,
                "local_path": final_path,
                "url":        pathway.diagram_url,
            })
    return items


async def _download_kegg_images(
    record,
    save_dir: str,
    annotate: bool,
    max_callouts: int,
    no_ai_polish: bool,
    annotate_dir: str,
) -> list[dict]:
    from fetchers.kegg import download_pathway_image
    ncbi_gene_id = getattr(record, "_ncbi_gene_id", "") or ""
    gene_symbols = [record.gene_name]
    items = []
    seen_pids: set[str] = set()
    for pathway in (getattr(record, "metabolic_pathways", []) +
                    getattr(record, "disease_pathways", [])):
        pid = getattr(pathway, "pathway_id", "")
        if not pid or pid in seen_pids:
            continue
        seen_pids.add(pid)
        plain_path, colored_path = await download_pathway_image(pid, ncbi_gene_id, save_dir)
        if colored_path:
            final_path = await _maybe_annotate(
                colored_path, pid, getattr(pathway, "pathway_name", pid),
                "KEGG", getattr(pathway, "diagram_url", ""), gene_symbols,
                record, annotate, max_callouts, no_ai_polish, annotate_dir,
            )
            pathway.diagram_image_path = final_path
            items.append({
                "pathway_id": pid,
                "local_path": final_path,
                "url":        getattr(pathway, "diagram_url", ""),
            })
    return items


async def _download_reactome_images(
    record,
    save_dir: str,
    annotate: bool,
    max_callouts: int,
    no_ai_polish: bool,
    annotate_dir: str,
) -> list[dict]:
    from fetchers.reactome import download_pathway_diagram
    uniprot_id = record.uniprot_id
    gene_symbols = [record.gene_name]
    items = []
    seen_pids: set[str] = set()
    for pathway in (getattr(record, "signaling_pathways", []) +
                    getattr(record, "disease_pathways", [])):
        pid    = getattr(pathway, "pathway_id", "")
        source = getattr(pathway, "pathway_source", "")
        if not pid or source != "Reactome" or pid in seen_pids:
            continue
        seen_pids.add(pid)
        plain_path, colored_path = await download_pathway_diagram(pid, uniprot_id, save_dir)
        if colored_path:
            final_path = await _maybe_annotate(
                colored_path, pid, getattr(pathway, "pathway_name", pid),
                "Reactome", getattr(pathway, "diagram_url", ""), gene_symbols,
                record, annotate, max_callouts, no_ai_polish, annotate_dir,
            )
            pathway.diagram_image_path = final_path
            items.append({
                "pathway_id": pid,
                "local_path": final_path,
                "url":        getattr(pathway, "diagram_url", ""),
            })
    return items


async def _download_string_image(
    record,
    save_dir: str,
    annotate: bool,
    max_callouts: int,
    no_ai_polish: bool,
    annotate_dir: str,
) -> list[dict]:
    from fetchers.string_db import download_network_image
    gene_symbol = record.gene_name
    local_path  = await download_network_image(gene_symbol, save_dir)
    if local_path:
        string_url = (
            f"https://string-db.org/cgi/network?"
            f"identifiers={gene_symbol}&species=9606"
        )
        final_path = await _maybe_annotate(
            local_path, f"string_{gene_symbol}",
            f"{gene_symbol} interaction network",
            "STRING", string_url, [gene_symbol],
            record, annotate, max_callouts, no_ai_polish, annotate_dir,
        )
        return [{"gene": gene_symbol, "local_path": final_path, "url": string_url}]
    return []
