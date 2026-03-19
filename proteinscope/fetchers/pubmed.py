"""PubMed / NCBI E-utilities fetcher for full citation metadata.

Given a list of PubMed IDs, fetches title, authors, journal, year, and DOI
via the Efetch XML endpoint and parses with lxml.
"""

from __future__ import annotations

import os
from typing import Optional

import httpx
from lxml import etree

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _api_key_params() -> dict:
    key = os.getenv("NCBI_API_KEY")
    return {"api_key": key} if key else {}


async def fetch_pubmed_xml(pmids: list[str], client: httpx.AsyncClient) -> bytes:
    """Fetch PubMed records in XML for a list of PMIDs."""
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        **_api_key_params(),
    }
    r = await client.get(f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=60)
    r.raise_for_status()
    return r.content


def parse_pubmed_xml(xml_bytes: bytes) -> list[dict]:
    """Parse PubMed XML response into list of citation dicts."""
    root = etree.fromstring(xml_bytes)
    refs = []
    for article in root.findall(".//PubmedArticle"):
        medline = article.find("MedlineCitation")
        if medline is None:
            continue
        art = medline.find("Article")
        if art is None:
            continue

        # Title
        title_el = art.find("ArticleTitle")
        title = (title_el.text or "") if title_el is not None else ""

        # Authors
        authors = []
        for author in art.findall(".//Author"):
            last = author.findtext("LastName", "")
            initials = author.findtext("Initials", "")
            if last:
                authors.append(f"{last} {initials}".strip())

        # Journal
        journal_el = art.find("Journal")
        journal = ""
        year = 0
        if journal_el is not None:
            journal = journal_el.findtext("Title", "") or journal_el.findtext("ISOAbbreviation", "")
            pub_date = journal_el.find(".//PubDate")
            if pub_date is not None:
                year_text = pub_date.findtext("Year", "0")
                try:
                    year = int(year_text)
                except ValueError:
                    year = 0

        # PubMed ID
        pmid = medline.findtext("PMID", "")

        # DOI
        doi = None
        for id_el in article.findall(".//ArticleId"):
            if id_el.get("IdType") == "doi":
                doi = id_el.text
                break

        refs.append({
            "pubmed_id": pmid,
            "doi": doi,
            "title": title,
            "authors": authors,
            "journal": journal,
            "year": year,
        })
    return refs


async def fetch_citations(pmids: list[str]) -> list[dict]:
    """Top-level: given PMID list, return citation dicts (batched to avoid URI limits)."""
    if not pmids:
        return []
    async with httpx.AsyncClient() as client:
        results = []
        batch_size = 20
        for i in range(0, len(pmids), batch_size):
            batch = pmids[i: i + batch_size]
            try:
                xml_bytes = await fetch_pubmed_xml(batch, client)
                results.extend(parse_pubmed_xml(xml_bytes))
            except Exception:
                continue
        return results
