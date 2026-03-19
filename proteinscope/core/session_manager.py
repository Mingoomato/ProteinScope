"""Session Context Layer (Phase 4).

Enables pronoun resolution ("now compare it to mouse"), multi-protein
workflows, and context injection into Gemini prompts.

SQLite-backed storage in data/sessions.db (auto-created).
Each session holds the last 20 protein analyses; oldest entries are
dropped when the limit is exceeded.
"""

from __future__ import annotations

import json
import logging
import re
import sqlite3
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from pydantic import BaseModel, Field

_log = logging.getLogger(__name__)

# DB path — created alongside the proteinscope package root
_DB_PATH = Path(__file__).parent.parent / "data" / "sessions.db"
_DB_PATH.parent.mkdir(parents=True, exist_ok=True)
_MAX_HISTORY = 20


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------


class SessionHistoryEntry(BaseModel):
    job_id: str
    gene: str
    accession: str
    timestamp: str          # ISO 8601
    key_findings: str       # 1-sentence Gemini summary (or fallback)


class SessionContext(BaseModel):
    session_id: str
    created_at: str         # ISO 8601
    proteins_analyzed: List[str] = Field(default_factory=list)    # UniProt accessions in order
    named_entities: Dict[str, str] = Field(default_factory=dict)  # alias → accession
    last_job_id: Optional[str] = None
    history: List[SessionHistoryEntry] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------


def _get_conn() -> sqlite3.Connection:
    conn = sqlite3.connect(str(_DB_PATH), check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


def _ensure_schema(conn: sqlite3.Connection) -> None:
    conn.execute("""
        CREATE TABLE IF NOT EXISTS sessions (
            session_id TEXT PRIMARY KEY,
            data       TEXT NOT NULL
        )
    """)
    conn.commit()


# ---------------------------------------------------------------------------
# Public SessionManager class
# ---------------------------------------------------------------------------


class SessionManager:
    """Thin wrapper around an SQLite sessions table."""

    def get_or_create(self, session_id: str) -> SessionContext:
        """Return existing SessionContext or create a new one."""
        try:
            conn = _get_conn()
            _ensure_schema(conn)
            row = conn.execute(
                "SELECT data FROM sessions WHERE session_id = ?", (session_id,)
            ).fetchone()
            conn.close()
            if row:
                return SessionContext(**json.loads(row["data"]))
        except Exception as exc:
            _log.debug("Session load failed for %s: %s", session_id, exc)

        ctx = SessionContext(
            session_id=session_id,
            created_at=datetime.utcnow().isoformat(),
        )
        self._save(ctx)
        return ctx

    def update(self, session_id: str, job_result: dict) -> None:
        """Update session with a completed job result."""
        try:
            ctx = self.get_or_create(session_id)

            gene = job_result.get("gene_name", "")
            accession = job_result.get("uniprot_id", "")
            job_id = job_result.get("job_id", "")
            organism = job_result.get("organism", "")

            if not gene or not accession:
                return

            # Build a brief key_findings string (use ai_summary if available)
            ai_sum = job_result.get("ai_summary") or ""
            if ai_sum:
                # First sentence only
                key_findings = ai_sum.split(".")[0].strip()[:200]
            else:
                protein_name = job_result.get("protein_name", gene)
                key_findings = f"{gene} ({protein_name}) in {organism}"

            entry = SessionHistoryEntry(
                job_id=job_id,
                gene=gene,
                accession=accession,
                timestamp=datetime.utcnow().isoformat(),
                key_findings=key_findings,
            )

            ctx.history.append(entry)
            if len(ctx.history) > _MAX_HISTORY:
                ctx.history = ctx.history[-_MAX_HISTORY:]

            if accession not in ctx.proteins_analyzed:
                ctx.proteins_analyzed.append(accession)

            ctx.last_job_id = job_id

            # Update named_entities
            ctx.named_entities["it"] = accession
            ctx.named_entities["the protein"] = accession
            ctx.named_entities[gene.lower()] = accession
            ctx.named_entities[protein_name.lower() if (protein_name := job_result.get("protein_name", "")) else gene.lower()] = accession

            self._save(ctx)

        except Exception as exc:
            _log.debug("Session update failed for %s: %s", session_id, exc)

    def get_gemini_context_block(self, session_id: str) -> str:
        """Return a context block string to prepend to Gemini prompts."""
        try:
            ctx = self.get_or_create(session_id)
            if not ctx.history:
                return ""

            lines = ["<session_context>"]
            lines.append("User has analyzed in this session:")
            for entry in ctx.history[-5:]:  # last 5 only to save tokens
                lines.append(f"  - {entry.gene} ({entry.accession}): {entry.key_findings}")

            if ctx.named_entities:
                relevant = {k: v for k, v in ctx.named_entities.items()
                            if k not in ("it", "the protein")}
                if relevant:
                    lines.append(f"Named proteins: {', '.join(f'{k}={v}' for k, v in list(relevant.items())[:5])}")

            lines.append("</session_context>")
            return "\n".join(lines)
        except Exception:
            return ""

    def resolve_entity(self, session_id: str, text: str) -> str:
        """Replace pronouns/aliases with canonical UniProt accessions in text."""
        try:
            ctx = self.get_or_create(session_id)
            if not ctx.named_entities:
                return text

            # Sort by length desc so longer aliases match first
            for alias, accession in sorted(
                ctx.named_entities.items(), key=lambda kv: -len(kv[0])
            ):
                if not alias:
                    continue
                pattern = re.compile(r"\b" + re.escape(alias) + r"\b", re.IGNORECASE)
                text = pattern.sub(accession, text, count=1)
            return text
        except Exception:
            return text

    def _save(self, ctx: SessionContext) -> None:
        try:
            conn = _get_conn()
            _ensure_schema(conn)
            conn.execute(
                "INSERT OR REPLACE INTO sessions (session_id, data) VALUES (?, ?)",
                (ctx.session_id, ctx.model_dump_json()),
            )
            conn.commit()
            conn.close()
        except Exception as exc:
            _log.debug("Session save failed: %s", exc)


# Module-level singleton
session_manager = SessionManager()
