"""Multimodal Embedding Search — FAISS concept search over cached proteins.

Uses sentence-transformers for text encoding and faiss-cpu for ANN search.
On first call, builds an index from cached protein data. Index stored at
data/protein_index.faiss; metadata at data/protein_index_meta.json.

Enables natural language concept queries like:
  "kinase with disordered transactivation domain"
  "membrane protein not expressed in liver"
"""

from __future__ import annotations
import json
import logging
import os
from pathlib import Path
from typing import List, Optional
import numpy as np

_log = logging.getLogger(__name__)
_INDEX_PATH = Path(__file__).parent.parent / "data" / "protein_index.faiss"
_META_PATH = Path(__file__).parent.parent / "data" / "protein_index_meta.json"
_INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)

_encoder = None  # lazy-loaded SentenceTransformer
_faiss_index = None
_meta: list[dict] = []   # [{gene, accession, organism, description}]


def _get_encoder():
    global _encoder
    if _encoder is None:
        from sentence_transformers import SentenceTransformer
        _encoder = SentenceTransformer("all-MiniLM-L6-v2")
    return _encoder


def _get_index_and_meta() -> tuple:
    """Load or return cached FAISS index + metadata."""
    global _faiss_index, _meta
    if _faiss_index is not None:
        return _faiss_index, _meta
    if _INDEX_PATH.exists() and _META_PATH.exists():
        try:
            import faiss
            _faiss_index = faiss.read_index(str(_INDEX_PATH))
            _meta = json.loads(_META_PATH.read_text(encoding="utf-8"))
            return _faiss_index, _meta
        except Exception as exc:
            _log.debug("Failed to load FAISS index: %s", exc)
    return None, []


def index_protein(gene: str, accession: str, organism: str, description: str,
                  go_terms: list[str] = None) -> None:
    """Add a single protein to the FAISS index (called after each successful search)."""
    global _faiss_index, _meta
    try:
        import faiss
        encoder = _get_encoder()
        text = f"{gene} {accession} {organism} {description} {' '.join(go_terms or [])}"
        vec = encoder.encode([text], normalize_embeddings=True).astype("float32")
        dim = vec.shape[1]

        # Check for duplicates
        if any(m["accession"] == accession for m in _meta):
            return

        if _faiss_index is None:
            _faiss_index = faiss.IndexFlatIP(dim)  # inner product = cosine on normalized vecs

        _faiss_index.add(vec)
        _meta.append({"gene": gene, "accession": accession,
                       "organism": organism, "description": description[:200]})

        # Persist
        faiss.write_index(_faiss_index, str(_INDEX_PATH))
        _META_PATH.write_text(json.dumps(_meta), encoding="utf-8")
    except Exception as exc:
        _log.debug("index_protein failed: %s", exc)


def search_by_concept(query: str, top_k: int = 10) -> list[dict]:
    """
    Search cached proteins by natural language concept.
    Returns ranked list of [{gene, accession, organism, description, score}].
    Returns [] if index is empty or packages not installed.
    """
    try:
        idx, meta = _get_index_and_meta()
        if idx is None or idx.ntotal == 0:
            return []
        encoder = _get_encoder()
        vec = encoder.encode([query], normalize_embeddings=True).astype("float32")
        k = min(top_k, idx.ntotal)
        scores, indices = idx.search(vec, k)
        results = []
        for score, i in zip(scores[0], indices[0]):
            if i < 0 or i >= len(meta):
                continue
            results.append({**meta[i], "score": float(score)})
        return results
    except Exception as exc:
        _log.debug("search_by_concept failed: %s", exc)
        return []


def get_index_size() -> int:
    """Return number of proteins currently indexed."""
    idx, _ = _get_index_and_meta()
    return idx.ntotal if idx else 0
