"""UniProtKB REST API v2 fetcher.

Retrieves protein identity, function, cofactors, subcellular location,
features (binding/active sites, domains), isoforms, and references.
"""

from __future__ import annotations

import httpx
from typing import Optional

BASE = "https://rest.uniprot.org/uniprotkb"


async def fetch_by_accession(accession: str, client: Optional[httpx.AsyncClient] = None) -> dict:
    """Fetch a full UniProt entry by accession ID."""
    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=30)
    try:
        r = await client.get(f"{BASE}/{accession}", params={"format": "json"})
        r.raise_for_status()
        return r.json()
    finally:
        if own_client:
            await client.aclose()


async def fetch_fasta(accession: str, client: Optional[httpx.AsyncClient] = None) -> str:
    """Fetch canonical FASTA sequence."""
    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=30)
    try:
        r = await client.get(f"{BASE}/{accession}.fasta")
        r.raise_for_status()
        return r.text
    finally:
        if own_client:
            await client.aclose()


async def search_protein(
    query: str,
    organism: Optional[str] = None,
    size: int = 5,
    client: Optional[httpx.AsyncClient] = None,
) -> list[dict]:
    """Search UniProt by protein name / gene symbol. Returns top results."""
    q = query
    if organism:
        q += f" AND organism_name:{organism}"
    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=30)
    try:
        r = await client.get(
            f"{BASE}/search",
            params={"query": q, "format": "json", "size": size},
        )
        r.raise_for_status()
        return r.json().get("results", [])
    finally:
        if own_client:
            await client.aclose()


def extract_protein_name(entry: dict) -> str:
    desc = entry.get("proteinDescription", {})
    rec = desc.get("recommendedName", {})
    full = rec.get("fullName", {})
    return full.get("value", "Unknown")


def extract_gene_name(entry: dict) -> str:
    genes = entry.get("genes", [])
    if genes:
        gn = genes[0].get("geneName", {})
        return gn.get("value", "")
    return ""


def extract_organism(entry: dict) -> tuple[str, str]:
    """Return (scientific_name, taxonomy_id)."""
    org = entry.get("organism", {})
    name = org.get("scientificName", "")
    tax_id = str(org.get("taxonId", ""))
    return name, tax_id


def extract_function(entry: dict) -> str:
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                return texts[0].get("value", "")
    return ""


def extract_cofactors(entry: dict) -> list[str]:
    cofactors = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "COFACTOR":
            for cf in comment.get("cofactors", []):
                name = cf.get("name", "")
                if name:
                    cofactors.append(name)
    return cofactors


def extract_subcellular_location(entry: dict) -> list[str]:
    locations = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for loc_block in comment.get("subcellularLocations", []):
                loc = loc_block.get("location", {})
                val = loc.get("value", "")
                if val:
                    locations.append(val)
    return locations


def extract_go_terms(entry: dict) -> tuple[list[str], list[str]]:
    """Return (biological_process_terms, molecular_function_terms)."""
    bp, mf = [], []
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GO":
            for prop in ref.get("properties", []):
                if prop.get("key") == "GoTerm":
                    val = prop.get("value", "")
                    if val.startswith("P:"):
                        bp.append(val[2:])
                    elif val.startswith("F:"):
                        mf.append(val[2:])
    return bp, mf


def extract_features(entry: dict, canonical_sequence: str) -> tuple[list[dict], list[dict], list[dict]]:
    """Return (binding_sites, active_sites, domains) as raw dicts with sequence fragments."""
    binding, active, domains = [], [], []
    for feat in entry.get("features", []):
        ftype = feat.get("type", "").lower()
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        # UniProt is 1-based inclusive
        fragment = canonical_sequence[start - 1: end] if canonical_sequence else ""
        desc = feat.get("description", "")
        ligand_raw = feat.get("ligand")
        if isinstance(ligand_raw, dict):
            ligand_name = ligand_raw.get("name", "") or ""
        elif isinstance(ligand_raw, str):
            ligand_name = ligand_raw
        else:
            # Fallback: try deprecated `ligands` list format
            ligand_name = ""
            for lig in feat.get("ligands", []):
                ligand_name = lig.get("name", "") if isinstance(lig, dict) else str(lig)
                break

        record = {
            "feature_type": feat.get("type", ""),
            "start": start,
            "end": end,
            "sequence_fragment": fragment,
            "description": desc or None,
            "ligand": ligand_name or None,
        }

        if ftype in ("binding site",):
            binding.append(record)
        elif ftype in ("active site",):
            active.append(record)
        elif ftype in ("domain", "region", "motif", "zinc finger"):
            domains.append(record)
    return binding, active, domains


def extract_references(entry: dict) -> list[dict]:
    """Extract PubMed IDs and basic metadata from UniProt references."""
    refs = []
    for ref in entry.get("references", []):
        citation = ref.get("citation", {})
        pubmed_id = None
        doi = None
        for db_ref in citation.get("citationCrossReferences", []):
            if db_ref.get("database") == "PubMed":
                pubmed_id = db_ref.get("id")
            elif db_ref.get("database") == "DOI":
                doi = db_ref.get("id")
        if pubmed_id or doi:
            authors = [
                a if isinstance(a, str) else a.get("value", "")
                for a in citation.get("authors", [])
            ]
            refs.append({
                "pubmed_id": pubmed_id,
                "doi": doi,
                "title": citation.get("title", ""),
                "authors": authors,
                "journal": citation.get("journal", ""),
                "year": citation.get("publicationDate", "0")[:4],
            })
    return refs


def extract_isoforms_from_comments(entry: dict) -> list[dict]:
    """Extract isoform metadata from ALTERNATIVE_SEQUENCE comments."""
    isoforms = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "ALTERNATIVE SEQUENCE":
            for iso in comment.get("isoforms", []):
                iso_id = iso.get("id", {})
                if isinstance(iso_id, list):
                    iso_id = iso_id[0] if iso_id else ""
                name_block = iso.get("name", {})
                name = name_block.get("value", iso_id) if isinstance(name_block, dict) else iso_id
                note_block = iso.get("note", {})
                note = note_block.get("value", "") if isinstance(note_block, dict) else ""
                isoforms.append({
                    "isoform_id": iso_id,
                    "name": name,
                    "differences": note,
                })
    return isoforms


def get_ncbi_gene_id(entry: dict) -> Optional[str]:
    """Pull the NCBI GeneID cross-reference from UniProt entry."""
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GeneID":
            return ref.get("id")
    return None
