"""Isoform handler.

Parses UniProt ALTERNATIVE_SEQUENCE features and isoform FASTA entries.
Reconstructs each isoform sequence and generates a diff description
against the canonical.
"""

from __future__ import annotations

import re
from typing import Optional

import httpx

from core.models import Isoform

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


async def fetch_isoform_fastas(accession: str) -> dict[str, str]:
    """Fetch all isoform FASTA sequences for an accession.

    Returns dict of {isoform_id: sequence}.
    """
    async with httpx.AsyncClient(timeout=30) as client:
        r = await client.get(
            f"{UNIPROT_BASE}/{accession}.fasta",
            params={"includeIsoform": "true"},
        )
        if r.status_code != 200:
            return {}
        return _parse_fasta_block(r.text)


def _parse_fasta_block(fasta_text: str) -> dict[str, str]:
    """Parse a multi-FASTA string into {id: sequence} dict."""
    seqs: dict[str, str] = {}
    current_id: Optional[str] = None
    current_seq: list[str] = []

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if current_id:
                seqs[current_id] = "".join(current_seq)
            # Extract isoform accession from header, e.g. >sp|P68871-2|...
            header_id = line[1:].split()[0]
            current_id = header_id
            current_seq = []
        else:
            current_seq.append(line.strip())

    if current_id:
        seqs[current_id] = "".join(current_seq)

    return seqs


def build_isoforms(
    uniprot_entry: dict,
    canonical_sequence: str,
    isoform_seqs: dict[str, str],
) -> list[Isoform]:
    """Build Isoform objects from UniProt entry and fetched FASTA sequences."""
    isoforms: list[Isoform] = []

    # Collect isoform metadata from ALTERNATIVE SEQUENCE comments
    iso_meta: dict[str, dict] = {}
    for comment in uniprot_entry.get("comments", []):
        if comment.get("commentType") != "ALTERNATIVE SEQUENCE":
            continue
        for iso in comment.get("isoforms", []):
            iso_ids = iso.get("id", [])
            if isinstance(iso_ids, str):
                iso_ids = [iso_ids]
            for iso_id in iso_ids:
                name_block = iso.get("name", {})
                name = name_block.get("value", iso_id) if isinstance(name_block, dict) else str(iso_id)
                note_block = iso.get("note", {})
                note = note_block.get("value", "") if isinstance(note_block, dict) else ""
                iso_meta[iso_id] = {"name": name, "note": note}

    # Build Isoform objects
    accession = uniprot_entry.get("primaryAccession", "")
    for iso_id, seq in isoform_seqs.items():
        # Skip canonical (no dash in ID or matches base accession)
        if "-" not in iso_id and iso_id == accession:
            continue
        meta = iso_meta.get(iso_id, {})
        name = meta.get("name", iso_id)
        note = meta.get("note", "")
        if not note:
            note = _simple_diff(canonical_sequence, seq)
        isoforms.append(Isoform(
            isoform_id=iso_id,
            name=name,
            sequence=seq,
            differences=note,
        ))

    return isoforms


def _simple_diff(canonical: str, isoform: str) -> str:
    """Generate a simple length-based diff description."""
    len_diff = len(isoform) - len(canonical)
    if len_diff == 0:
        return "Same length as canonical; internal substitution(s) present."
    elif len_diff > 0:
        return f"{len_diff} residues longer than canonical."
    else:
        return f"{abs(len_diff)} residues shorter than canonical."
