"""UniProt proteome fetcher for Target Discovery Mode.

Fetches reviewed (Swiss-Prot) UniProt entries for a given organism name.
Used by analyzers/target_discovery_analyzer.py to build the candidate protein list.

Returns empty list on any failure — graceful degradation guaranteed.
"""
from __future__ import annotations

import logging
from typing import Optional
import httpx

_log = logging.getLogger(__name__)
_BASE = "https://rest.uniprot.org/uniprotkb/search"


async def fetch_organism_proteome(
    organism_name: str,
    max_proteins: int = 80,
) -> list[dict]:
    """Fetch UniProt entries for a given organism name.

    Tries reviewed (Swiss-Prot) entries first; falls back to including
    unreviewed (TrEMBL) entries if fewer than 5 results are returned.

    Args:
        organism_name: Scientific or common name, e.g. "Malassezia globosa"
        max_proteins:  Maximum number of entries to return (default 80)

    Returns:
        List of dicts with keys:
            accession (str)
            gene_name (str)
            protein_name (str)
            function_summary (str)
            go_terms (list[str])
            sequence_length (int)
    """
    fields = "accession,gene_names,protein_name,cc_function,go,length"

    async def _fetch(query: str) -> list[dict]:
        try:
            async with httpx.AsyncClient(timeout=20.0) as client:
                resp = await client.get(
                    _BASE,
                    params={
                        "query": query,
                        "fields": fields,
                        "format": "json",
                        "size": max_proteins,
                    },
                )
                if resp.status_code != 200:
                    return []
                data = resp.json()
                results = []
                for entry in (data.get("results") or []):
                    acc = entry.get("primaryAccession", "")
                    # Gene name
                    gene_names_block = entry.get("genes") or []
                    gene_name = ""
                    if gene_names_block:
                        gene_name = (
                            gene_names_block[0].get("geneName", {}).get("value", "")
                            or (gene_names_block[0].get("synonyms") or [{}])[0].get("value", "")
                        )
                    # Protein name
                    pn_block = entry.get("proteinDescription", {})
                    protein_name = (
                        pn_block.get("recommendedName", {}).get("fullName", {}).get("value", "")
                        or pn_block.get("submittedName", [{}])[0].get("fullName", {}).get("value", "")
                        if pn_block
                        else ""
                    )
                    # Function summary from comments
                    function_summary = ""
                    for comment in (entry.get("comments") or []):
                        if comment.get("commentType") == "FUNCTION":
                            texts = comment.get("texts") or []
                            if texts:
                                function_summary = texts[0].get("value", "")
                                break
                    # GO terms
                    go_terms: list[str] = []
                    for ref in (entry.get("uniProtKBCrossReferences") or []):
                        if ref.get("database") == "GO":
                            for prop in (ref.get("properties") or []):
                                if prop.get("key") == "GoTerm":
                                    go_terms.append(prop.get("value", ""))
                    # Sequence length
                    seq_length = entry.get("sequence", {}).get("length", 0)

                    results.append({
                        "accession": acc,
                        "gene_name": gene_name,
                        "protein_name": protein_name,
                        "function_summary": function_summary,
                        "go_terms": go_terms[:10],  # cap at 10
                        "sequence_length": seq_length,
                    })
                return results
        except Exception as exc:
            _log.debug("fetch_organism_proteome error for %s: %s", organism_name, exc)
            return []

    # Try reviewed entries first
    organism_escaped = organism_name.replace('"', '')
    query_reviewed = f'organism_name:"{organism_escaped}" AND reviewed:true'
    results = await _fetch(query_reviewed)

    # Fallback: include unreviewed if fewer than 5 results
    if len(results) < 5:
        query_all = f'organism_name:"{organism_escaped}"'
        results = await _fetch(query_all)
        _log.debug(
            "Proteome fallback (unreviewed included) for %s: %d results",
            organism_name,
            len(results),
        )

    return results
