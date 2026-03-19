"""Literature RAG — Retrieval-Augmented Generation over PubMed abstracts.

Uses sentence-transformers (all-MiniLM-L6-v2) to embed abstracts and
faiss-cpu for vector similarity search.

Workflow:
1. LiteratureRAG.index(references) — embeds abstracts into FAISS index
2. LiteratureRAG.query(question, top_k) — retrieves most relevant abstracts
3. LiteratureRAG.answer(question) — calls Gemini with retrieved context

The model is downloaded once and cached in ~/.cache/huggingface/.
"""

from __future__ import annotations

import logging

_log = logging.getLogger(__name__)


class LiteratureRAG:
    """Retrieval-Augmented Generation engine over a set of protein literature references."""

    def __init__(self) -> None:
        self._model = None   # lazy-loaded SentenceTransformer
        self._index = None   # faiss.IndexFlatIP
        self._refs: list = []  # stored references
        self._device: str | None = None  # resolved on first model load

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _get_model(self):
        """Lazy-load all-MiniLM-L6-v2 on GPU if available, otherwise CPU.

        Uses utils.device.get_device() so the GPU selection policy is shared
        with all other PyTorch-dependent modules in this package.
        """
        if self._model is None:
            try:
                from sentence_transformers import SentenceTransformer
                from utils.device import get_device
                self._device = get_device()
                _log.info("Loading SentenceTransformer on device=%s", self._device)
                self._model = SentenceTransformer(
                    "all-MiniLM-L6-v2",
                    device=self._device,
                )
            except Exception:
                pass
        return self._model

    def _embed(self, texts: list[str]):
        """Embed a list of texts. Returns a numpy float32 array or None."""
        model = self._get_model()
        if model is None or not texts:
            return None
        try:
            import numpy as np
            embeddings = model.encode(texts, show_progress_bar=False, convert_to_numpy=True)
            return embeddings.astype("float32")
        except Exception:
            return None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def index(self, references: list) -> int:
        """Embed and index reference abstracts.

        Args:
            references: List of core.models.Reference objects (or any objects with
                        title, authors, journal, year, apa_citation attributes).

        Returns:
            Number of successfully indexed references (0 if model unavailable).
        """
        if not references:
            return 0

        # Build text corpus: concatenate key bibliographic fields
        texts: list[str] = []
        for ref in references:
            try:
                authors_str = " ".join(ref.authors[:3]) if ref.authors else ""
                text = f"{ref.title} {authors_str} {ref.journal} {ref.year}"
                texts.append(text.strip())
            except Exception:
                texts.append("")

        vectors = self._embed(texts)
        if vectors is None:
            return 0

        try:
            import faiss
            import numpy as np

            faiss.normalize_L2(vectors)
            index = faiss.IndexFlatIP(vectors.shape[1])
            index.add(vectors)

            self._index = index
            self._refs = list(references)
            return len(references)
        except Exception:
            return 0

    def query(self, question: str, top_k: int = 3) -> list[dict]:
        """Find the most relevant references for a question.

        Args:
            question: Natural language question about the protein.
            top_k:    Maximum number of references to return.

        Returns:
            List of dicts with keys: reference, score, abstract_snippet.
            Returns empty list if index is not built or model unavailable.
        """
        if self._index is None or not self._refs:
            return []

        vectors = self._embed([question])
        if vectors is None:
            return []

        try:
            import faiss
            import numpy as np

            query_vec = vectors  # already float32 from _embed
            faiss.normalize_L2(query_vec)

            k = min(top_k, len(self._refs))
            distances, indices = self._index.search(query_vec, k)

            results: list[dict] = []
            for score, idx in zip(distances[0], indices[0]):
                if idx < 0 or idx >= len(self._refs):
                    continue
                ref = self._refs[idx]
                # Build a short abstract snippet from available metadata
                try:
                    snippet = f"{ref.title} ({ref.year}). {ref.journal}."
                except Exception:
                    snippet = str(ref)
                results.append({
                    "reference": ref,
                    "score": float(score),
                    "abstract_snippet": snippet,
                })
            return results
        except Exception:
            return []

    async def answer(self, question: str, top_k: int = 3) -> str:
        """Answer a question about the indexed protein using RAG + Gemini.

        Retrieves the most relevant references, builds a context block, and
        calls Gemini 2.5 Pro to answer the question with in-text citations.

        Args:
            question: The user's research question.
            top_k:    Number of references to include as context.

        Returns:
            Answer string with citations, or empty string if model/key unavailable.
        """
        hits = self.query(question, top_k=top_k)
        if not hits:
            # No index or no results — still attempt Gemini if key available,
            # but return empty string gracefully if Gemini also unavailable.
            context_block = "No literature context available."
        else:
            lines = []
            for i, hit in enumerate(hits, start=1):
                ref = hit["reference"]
                try:
                    citation = ref.apa_citation or f"{ref.title} ({ref.year}). {ref.journal}."
                except Exception:
                    citation = str(ref)
                lines.append(f"Paper {i}: {citation}")
            context_block = "\n".join(lines)

        prompt = (
            "You are a molecular biology expert answering a research question based on "
            "provided literature references.\n\n"
            "Literature context:\n"
            f"{context_block}\n\n"
            f"Question: {question}\n\n"
            "Answer the question concisely, citing papers by number (e.g. [1], [2]). "
            "If the context does not contain enough information, say so clearly. "
            "Do not invent data not present in the references above."
        )

        try:
            from core.gemini_interpreter import _call
            result = await _call(prompt)
            return result or ""
        except Exception:
            return ""
