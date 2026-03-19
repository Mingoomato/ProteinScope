"""NCBI E-utilities fetcher for DNA/CDS sequences.

Workflow:
1. Use NCBI Gene ID (from UniProt cross-reference) to find RefSeq mRNA via Elink.
2. Fetch GenBank record via Efetch.
3. Parse CDS feature with Biopython to extract nucleotide sequence.
4. Detect spliced exons (join(...) location syntax).
"""

from __future__ import annotations

import os
from io import StringIO
from typing import Optional

import httpx
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _api_key_params() -> dict:
    key = os.getenv("NCBI_API_KEY")
    return {"api_key": key} if key else {}


async def gene_id_to_mrna_ids(gene_id: str, client: httpx.AsyncClient) -> list[str]:
    """Use Elink to go from Gene ID -> Nucleotide (RefSeq mRNA) IDs."""
    params = {
        "dbfrom": "gene",
        "db": "nuccore",
        "id": gene_id,
        "linkname": "gene_nuccore_refseqrna",
        "retmode": "json",
        **_api_key_params(),
    }
    r = await client.get(f"{EUTILS_BASE}/elink.fcgi", params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    ids = []
    try:
        link_sets = data["linksets"][0]["linksetdbs"]
        for ls in link_sets:
            ids.extend(ls.get("links", []))
    except (KeyError, IndexError):
        pass
    return [str(i) for i in ids[:5]]  # limit to first 5


async def fetch_genbank_record(nuccore_id: str, client: httpx.AsyncClient) -> str:
    """Fetch GenBank text for a nucleotide record."""
    params = {
        "db": "nuccore",
        "id": nuccore_id,
        "rettype": "gb",
        "retmode": "text",
        **_api_key_params(),
    }
    r = await client.get(f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=60)
    r.raise_for_status()
    return r.text


def parse_cds_from_genbank(gb_text: str) -> tuple[Optional[str], bool]:
    """Parse CDS nucleotide sequence from GenBank text.

    Returns (cds_sequence, is_spliced).
    """
    handle = StringIO(gb_text)
    for record in SeqIO.parse(handle, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                is_spliced = isinstance(feature.location, CompoundLocation)
                cds_seq = str(feature.extract(record.seq))
                return cds_seq, is_spliced
    return None, False


async def fetch_cds_for_gene(gene_id: str) -> tuple[Optional[str], bool]:
    """Top-level: given an NCBI Gene ID, return (cds_sequence, is_spliced)."""
    async with httpx.AsyncClient() as client:
        mrna_ids = await gene_id_to_mrna_ids(gene_id, client)
        if not mrna_ids:
            return None, False
        gb_text = await fetch_genbank_record(mrna_ids[0], client)
        return parse_cds_from_genbank(gb_text)
