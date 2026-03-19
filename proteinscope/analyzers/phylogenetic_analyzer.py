"""Phylogenetic analyzer — build NJ tree from ortholog sequences.

Uses BioPython's Phylo module and MAFFT for MSA (falls back to pairwise).
Fetches ortholog sequences from UniProt search.
"""

from __future__ import annotations

import io
import math
from typing import Optional

import httpx
from pydantic import BaseModel, Field

from analyzers.sequence_aligner import build_msa, align_pairwise

# ---------------------------------------------------------------------------
# UniProt REST API endpoints
# ---------------------------------------------------------------------------
_UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class PhyloResult(BaseModel):
    newick_tree: str
    species_list: list[str]
    n_sequences: int
    evolutionary_rates: list[float]  # per position; low = conserved
    mean_identity_pct: float


# ---------------------------------------------------------------------------
# UniProt ortholog fetcher
# ---------------------------------------------------------------------------

async def fetch_orthologs(
    gene_name: str,
    n_species: int = 10,
) -> dict[str, str]:
    """Fetch reviewed UniProt entries for a gene across multiple organisms.

    Searches UniProt for ``gene_exact:{gene_name} AND reviewed:true`` and
    returns up to ``n_species`` results as a mapping of organism → sequence.

    Args:
        gene_name:  HGNC gene symbol (e.g. ``TP53``, ``EGFR``).
        n_species:  Maximum number of species sequences to return.

    Returns:
        Dict of {organism_name: amino_acid_sequence}. Empty dict on failure.
    """
    if not gene_name:
        return {}

    try:
        params = {
            "query": f"gene_exact:{gene_name} AND reviewed:true",
            "format": "fasta",
            "size": str(n_species),
            "fields": "sequence,organism_name",
        }

        async with httpx.AsyncClient(timeout=30) as client:
            response = await client.get(_UNIPROT_SEARCH, params=params)
            if response.status_code != 200:
                return {}

            fasta_text = response.text

        # Parse FASTA — use organism name from header as key
        result: dict[str, str] = {}
        current_header: Optional[str] = None
        buf: list[str] = []

        for line in fasta_text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header and buf:
                    result[current_header] = "".join(buf)
                # Extract organism name from header: "OS=Homo sapiens"
                header = line[1:]
                organism = _extract_organism(header)
                accession = header.split()[0] if header else "unknown"
                # Use organism+accession as key to handle duplicates
                current_header = f"{organism} ({accession})" if organism else accession
                buf = []
            else:
                buf.append(line)

        if current_header and buf:
            result[current_header] = "".join(buf)

        return result

    except Exception:
        return {}


def _extract_organism(fasta_header: str) -> str:
    """Extract OS= organism name from a UniProt FASTA header."""
    try:
        os_start = fasta_header.find("OS=")
        if os_start == -1:
            return ""
        os_end = fasta_header.find(" OX=", os_start)
        if os_end == -1:
            os_end = fasta_header.find(" GN=", os_start)
        if os_end == -1:
            os_end = len(fasta_header)
        return fasta_header[os_start + 3:os_end].strip()
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# Tree construction
# ---------------------------------------------------------------------------

def build_tree_from_sequences(sequences: dict[str, str]) -> Optional[PhyloResult]:
    """Construct a neighbor-joining phylogenetic tree from protein sequences.

    Steps:
    1. Build MSA (MAFFT → BioPython pairwise fallback).
    2. Compute pairwise identity matrix.
    3. Construct NJ tree with BioPython DistanceTreeConstructor.
    4. Extract Newick string and per-column evolutionary rates.

    Args:
        sequences: Mapping of {label: sequence}.

    Returns:
        PhyloResult with Newick string, species list, and evolutionary rates,
        or None if fewer than 2 sequences are provided or construction fails.
    """
    if not sequences or len(sequences) < 2:
        return None

    try:
        from Bio.Phylo.TreeConstruction import (
            DistanceMatrix,
            DistanceTreeConstructor,
        )
        from Bio import Phylo
        import numpy as np

        names = list(sequences.keys())
        seqs = list(sequences.values())
        n = len(names)

        # ── Step 1: MSA ───────────────────────────────────────────────────
        msa_result = build_msa(sequences)

        # ── Step 2: Pairwise identity matrix ──────────────────────────────
        # identity_matrix[i][j] = identity fraction (0–1) between i and j
        identity_matrix = _compute_identity_matrix(names, seqs)

        # Convert to distance matrix (lower-triangular for BioPython)
        # BioPython DistanceMatrix requires a lower-triangular list of lists
        matrix_data: list[list[float]] = []
        for i in range(n):
            row: list[float] = []
            for j in range(i + 1):
                if i == j:
                    row.append(0.0)
                else:
                    # distance = 1 - identity
                    row.append(1.0 - identity_matrix[i][j])
            matrix_data.append(row)

        dm = DistanceMatrix(names=names, matrix=matrix_data)

        # ── Step 3: NJ tree ───────────────────────────────────────────────
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

        # ── Step 4: Newick string ──────────────────────────────────────────
        newick_io = io.StringIO()
        Phylo.write(tree, newick_io, "newick")
        newick_str = newick_io.getvalue().strip()

        # ── Step 5: Evolutionary rates from MSA conservation ──────────────
        if msa_result is not None and msa_result.conservation_scores:
            # evolutionary rate = 1 - conservation (high conservation = low rate)
            evo_rates = [
                round(1.0 - score, 4)
                for score in msa_result.conservation_scores
            ]
        else:
            evo_rates = []

        # ── Step 6: Mean pairwise identity ────────────────────────────────
        pairwise_identities: list[float] = []
        for i in range(n):
            for j in range(i + 1, n):
                pairwise_identities.append(identity_matrix[i][j] * 100.0)

        mean_identity = (
            sum(pairwise_identities) / len(pairwise_identities)
            if pairwise_identities
            else 0.0
        )

        return PhyloResult(
            newick_tree=newick_str,
            species_list=names,
            n_sequences=n,
            evolutionary_rates=evo_rates,
            mean_identity_pct=round(mean_identity, 2),
        )

    except ImportError:
        return None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _compute_identity_matrix(names: list[str], seqs: list[str]) -> list[list[float]]:
    """Compute n×n pairwise identity matrix.

    Uses BioPython PairwiseAligner via align_pairwise() from sequence_aligner.
    Falls back to character-level identity for speed if alignment fails.
    """
    n = len(names)
    matrix: list[list[float]] = [[0.0] * n for _ in range(n)]

    for i in range(n):
        matrix[i][i] = 1.0
        for j in range(i + 1, n):
            result = align_pairwise(seqs[i], seqs[j])
            if result is not None:
                identity = result.identity_pct / 100.0
            else:
                # Naive fallback: character overlap
                min_len = min(len(seqs[i]), len(seqs[j]))
                if min_len == 0:
                    identity = 0.0
                else:
                    matches = sum(a == b for a, b in zip(seqs[i][:min_len], seqs[j][:min_len]))
                    identity = matches / max(len(seqs[i]), len(seqs[j]))

            matrix[i][j] = identity
            matrix[j][i] = identity

    return matrix
