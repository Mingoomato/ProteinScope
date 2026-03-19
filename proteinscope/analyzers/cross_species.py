"""Cross-species interaction domain analyzer.

Given a domain fragment sequence, this module:
1. POSTs to the UniProt BLAST API.
2. Polls for results with exponential backoff.
3. Filters hits by taxonomic group of interest.
4. Fetches binding site annotations for each ortholog.
5. Compares binding site sequences and flags compatibility.
"""

from __future__ import annotations

import asyncio
import time
from typing import Optional

import httpx

from core.models import CrossSpeciesEntry

UNIPROT_BLAST_SUBMIT = "https://www.uniprot.org/blast/run"
UNIPROT_BLAST_STATUS = "https://www.uniprot.org/blast/status"
UNIPROT_BLAST_RESULTS = "https://www.uniprot.org/blast/results"
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"

# Default identity threshold (%) for flagging compatibility
DEFAULT_THRESHOLD = 80.0


async def submit_blast(sequence: str, client: httpx.AsyncClient) -> Optional[str]:
    """Submit a BLAST job to UniProt. Returns job ID or None on failure."""
    data = {
        "sequence": sequence,
        "database": "uniprotkb_swissprot",
        "taxId": "",
        "threshold": "0.001",
        "matrix": "BLOSUM62",
        "filtering": "false",
        "gapped": "true",
        "hits": "50",
    }
    r = await client.post(UNIPROT_BLAST_SUBMIT, data=data, timeout=30)
    if r.status_code not in (200, 201, 303):
        return None
    # Response may include Location header with job ID
    location = r.headers.get("location", "")
    if location:
        return location.split("/")[-1]
    return r.text.strip() or None


async def poll_blast(job_id: str, client: httpx.AsyncClient, timeout_s: int = 120) -> bool:
    """Poll BLAST job until complete. Returns True on success."""
    backoff = 2
    elapsed = 0
    while elapsed < timeout_s:
        await asyncio.sleep(backoff)
        elapsed += backoff
        backoff = min(backoff * 2, 30)
        try:
            r = await client.get(f"{UNIPROT_BLAST_STATUS}/{job_id}", timeout=15)
            if r.status_code == 200:
                status = r.text.strip()
                if status in ("RUNNING", "PENDING", "QUEUED"):
                    continue
                return status == "FINISHED"
        except Exception:
            continue
    return False


async def get_blast_results(job_id: str, client: httpx.AsyncClient) -> list[dict]:
    """Retrieve BLAST results for a completed job."""
    r = await client.get(
        f"{UNIPROT_BLAST_RESULTS}/{job_id}",
        params={"format": "json"},
        timeout=30,
    )
    if r.status_code != 200:
        return []
    data = r.json()
    return data.get("hits", [])


async def fetch_uniprot_features(accession: str, client: httpx.AsyncClient) -> dict:
    """Fetch UniProt entry for binding site comparison."""
    r = await client.get(f"{UNIPROT_BASE}/{accession}", params={"format": "json"}, timeout=30)
    if r.status_code != 200:
        return {}
    return r.json()


def extract_binding_fragment(entry: dict, canonical: str) -> Optional[str]:
    """Get first binding site fragment from a UniProt entry."""
    for feat in entry.get("features", []):
        if feat.get("type", "").lower() == "binding site":
            loc = feat.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            if start and end and canonical:
                return canonical[start - 1: end]
    return None


def compute_identity(seq_a: str, seq_b: str) -> float:
    """Compute rough sequence identity % between two same-length fragments."""
    if not seq_a or not seq_b:
        return 0.0
    min_len = min(len(seq_a), len(seq_b))
    max_len = max(len(seq_a), len(seq_b))
    matches = sum(a == b for a, b in zip(seq_a[:min_len], seq_b[:min_len]))
    return round(100.0 * matches / max_len, 1)


async def analyze_cross_species(
    domain_fragment: str,
    reference_fragment: str,
    threshold: float = DEFAULT_THRESHOLD,
) -> list[CrossSpeciesEntry]:
    """Run full cross-species analysis for a domain fragment.

    Args:
        domain_fragment:    Interaction domain sequence from query protein.
        reference_fragment: Same fragment used as baseline for identity comparison.
        threshold:          Minimum identity % to flag as compatible.

    Returns:
        List of CrossSpeciesEntry objects.
    """
    entries: list[CrossSpeciesEntry] = []

    async with httpx.AsyncClient(timeout=60) as client:
        job_id = await submit_blast(domain_fragment, client)
        if not job_id:
            return entries

        success = await poll_blast(job_id, client)
        if not success:
            return entries

        hits = await get_blast_results(job_id, client)

        for hit in hits[:20]:
            accession = hit.get("accession", "")
            organism = hit.get("organism", {}).get("scientificName", "")
            identity_pct = float(hit.get("identity", 0))

            # Skip the query organism itself
            if identity_pct == 100.0:
                continue

            # Fetch binding site of ortholog for direct domain comparison
            ortholog_entry = await fetch_uniprot_features(accession, client)
            ortholog_seq = ortholog_entry.get("sequence", {}).get("value", "")
            ortholog_binding = extract_binding_fragment(ortholog_entry, ortholog_seq)

            if ortholog_binding and reference_fragment:
                domain_identity = compute_identity(reference_fragment, ortholog_binding)
            else:
                domain_identity = identity_pct  # fall back to overall identity

            compatible = domain_identity >= threshold

            entries.append(CrossSpeciesEntry(
                organism=organism,
                uniprot_id=accession,
                identity_pct=domain_identity,
                interaction_domain_seq=ortholog_binding,
                compatible=compatible,
            ))

    return entries
