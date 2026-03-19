"""MetaCyc / BioCyc fetcher — detailed metabolic enzyme roles and reactions.

BioCyc/MetaCyc provides richer enzyme role descriptions than KEGG.
Free academic API key available at https://biocyc.org/

Authentication: HTTP Basic auth using BIOCYC_USERNAME / BIOCYC_PASSWORD from .env
"""

from __future__ import annotations

import os

import httpx
from dotenv import load_dotenv

load_dotenv()

_BASE = "https://websvc.biocyc.org"
_USER = os.getenv("BIOCYC_USERNAME", "")
_PASS = os.getenv("BIOCYC_PASSWORD", "")
_AUTH = (_USER, _PASS) if _USER and _PASS else None


async def fetch_gene_detail(gene_symbol: str) -> str:
    """Return full XML for a human gene entry in BioCyc."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/getxml",
            params={"id": f"HUMAN:{gene_symbol}", "detail": "full"},
            auth=_AUTH,
            timeout=30,
        )
    return r.text if r.status_code == 200 else ""


async def fetch_enzyme_reactions(gene_symbol: str) -> str:
    """Return XML listing reactions catalyzed by this enzyme in BioCyc."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/getxml",
            params={"id": f"HUMAN:{gene_symbol}", "type": "reactions"},
            auth=_AUTH,
            timeout=30,
        )
    return r.text if r.status_code == 200 else ""


async def fetch_pathway_detail(pathway_id: str) -> str:
    """Return full BioCyc XML for a pathway entry by its ID."""
    async with httpx.AsyncClient() as client:
        r = await client.get(
            f"{_BASE}/getxml",
            params={"id": pathway_id, "detail": "full"},
            auth=_AUTH,
            timeout=30,
        )
    return r.text if r.status_code == 200 else ""


def _parse_ec_numbers(xml_text: str) -> list[str]:
    """Extract EC numbers from BioCyc XML (simple regex approach)."""
    import re
    return re.findall(r"EC-\d+\.\d+\.\d+\.\d+", xml_text)


def _parse_reaction_names(xml_text: str) -> list[str]:
    """Extract reaction common names from BioCyc XML."""
    import re
    return re.findall(r"<common-name[^>]*>([^<]+)</common-name>", xml_text)


async def fetch_all_metacyc(gene_symbol: str) -> dict:
    """Top-level wrapper returning MetaCyc data for a gene symbol."""
    if not _AUTH:
        return {
            "available": False,
            "reason": "BIOCYC_USERNAME/BIOCYC_PASSWORD not set in .env",
            "ec_numbers": [],
            "reaction_names": [],
        }

    gene_xml = await fetch_gene_detail(gene_symbol)
    rxn_xml = await fetch_enzyme_reactions(gene_symbol)

    ec_numbers = _parse_ec_numbers(gene_xml + rxn_xml)
    reaction_names = _parse_reaction_names(rxn_xml)

    return {
        "available": True,
        "ec_numbers": list(set(ec_numbers)),
        "reaction_names": list(set(reaction_names)),
        "gene_xml": gene_xml,
        "reactions_xml": rxn_xml,
    }
