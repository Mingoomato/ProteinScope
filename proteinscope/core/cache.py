"""Local SQLite cache for API responses.

Keyed by (accession, data_source) with a configurable TTL (default 7 days).
This avoids redundant network calls during iterative report generation.
"""

from __future__ import annotations

import json
import sqlite3
import time
from pathlib import Path
from typing import Optional

DEFAULT_DB_PATH = Path(__file__).parent.parent / ".cache" / "proteinscope.db"
DEFAULT_TTL_DAYS = 7


class Cache:
    def __init__(self, db_path: Path = DEFAULT_DB_PATH):
        db_path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(db_path), check_same_thread=False)
        self._create_table()

    def _create_table(self) -> None:
        self._conn.execute(
            """
            CREATE TABLE IF NOT EXISTS cache (
                key       TEXT PRIMARY KEY,
                data      TEXT NOT NULL,
                expires   REAL NOT NULL
            )
            """
        )
        self._conn.commit()

    def _key(self, accession: str, source: str) -> str:
        return f"{accession}::{source}"

    def get(self, accession: str, source: str) -> Optional[dict]:
        """Return cached data or None if missing / expired."""
        key = self._key(accession, source)
        row = self._conn.execute(
            "SELECT data, expires FROM cache WHERE key = ?", (key,)
        ).fetchone()
        if row is None:
            return None
        data_str, expires = row
        if time.time() > expires:
            self.invalidate(accession, source)
            return None
        return json.loads(data_str)

    def set(self, accession: str, source: str, data: dict, ttl_days: int = DEFAULT_TTL_DAYS) -> None:
        """Store data with an expiry timestamp."""
        key = self._key(accession, source)
        expires = time.time() + ttl_days * 86400
        self._conn.execute(
            "INSERT OR REPLACE INTO cache (key, data, expires) VALUES (?, ?, ?)",
            (key, json.dumps(data, ensure_ascii=False), expires),
        )
        self._conn.commit()

    def invalidate(self, accession: str, source: str) -> None:
        """Remove a specific cache entry."""
        key = self._key(accession, source)
        self._conn.execute("DELETE FROM cache WHERE key = ?", (key,))
        self._conn.commit()

    def clear_expired(self) -> int:
        """Delete all expired entries. Returns count removed."""
        cursor = self._conn.execute(
            "DELETE FROM cache WHERE expires <= ?", (time.time(),)
        )
        self._conn.commit()
        return cursor.rowcount

    def close(self) -> None:
        self._conn.close()

    def __enter__(self) -> "Cache":
        return self

    def __exit__(self, *_) -> None:
        self.close()
