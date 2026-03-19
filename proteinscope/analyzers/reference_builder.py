"""Reference builder: formats citations in APA, Vancouver, and BibTeX styles."""

from __future__ import annotations

import re
import unicodedata
from core.models import Reference


def _slugify(text: str) -> str:
    """Create a safe BibTeX key from text."""
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    text = re.sub(r"[^\w]", "_", text)
    return text[:40].strip("_").lower()


def _apa_authors(authors: list[str]) -> str:
    """Format author list in APA style."""
    if not authors:
        return "Unknown"
    formatted = []
    for author in authors:
        parts = author.split()
        if len(parts) >= 2:
            last = parts[0]
            initials = ". ".join(p[0] for p in parts[1:]) + "."
            formatted.append(f"{last}, {initials}")
        else:
            formatted.append(author)
    if len(formatted) == 1:
        return formatted[0]
    if len(formatted) <= 6:
        return ", ".join(formatted[:-1]) + f", & {formatted[-1]}"
    return ", ".join(formatted[:6]) + ", ... & " + formatted[-1]


def _vancouver_authors(authors: list[str]) -> str:
    """Format author list in Vancouver style (max 6 then et al.)."""
    if not authors:
        return "Unknown"
    formatted = []
    for author in authors[:6]:
        parts = author.split()
        if len(parts) >= 2:
            last = parts[0]
            initials = "".join(p[0] for p in parts[1:])
            formatted.append(f"{last} {initials}")
        else:
            formatted.append(author)
    result = ", ".join(formatted)
    if len(authors) > 6:
        result += ", et al"
    return result


def build_apa(ref: Reference) -> str:
    """Generate APA 7th edition citation string."""
    authors = _apa_authors(ref.authors)
    year = f"({ref.year})" if ref.year else "(n.d.)"
    title = ref.title or "Untitled"
    journal = f"*{ref.journal}*" if ref.journal else ""
    doi_part = f" https://doi.org/{ref.doi}" if ref.doi else ""
    return f"{authors} {year}. {title}. {journal}.{doi_part}"


def build_vancouver(ref: Reference) -> str:
    """Generate Vancouver citation string."""
    authors = _vancouver_authors(ref.authors)
    year = str(ref.year) if ref.year else "n.d."
    title = ref.title or "Untitled"
    journal = ref.journal or ""
    doi_part = f" doi:{ref.doi}" if ref.doi else ""
    return f"{authors}. {title}. {journal}. {year};{doi_part}"


def build_bibtex(ref: Reference) -> str:
    """Generate BibTeX entry."""
    first_author = ref.authors[0].split()[0] if ref.authors else "Unknown"
    key = f"{_slugify(first_author)}{ref.year}"
    authors_field = " and ".join(ref.authors) if ref.authors else "Unknown"
    lines = [
        f"@article{{{key},",
        f"  author  = {{{authors_field}}},",
        f"  title   = {{{ref.title or 'Untitled'}}},",
        f"  journal = {{{ref.journal or ''}}},",
        f"  year    = {{{ref.year}}},",
    ]
    if ref.doi:
        lines.append(f"  doi     = {{{ref.doi}}},")
    if ref.pubmed_id:
        lines.append(f"  pmid    = {{{ref.pubmed_id}}},")
    lines.append("}")
    return "\n".join(lines)


def populate_citations(refs: list[Reference]) -> list[Reference]:
    """Populate apa_citation, vancouver_citation, and bibtex fields for each Reference."""
    populated = []
    for ref in refs:
        populated.append(ref.model_copy(update={
            "apa_citation": build_apa(ref),
            "vancouver_citation": build_vancouver(ref),
            "bibtex": build_bibtex(ref),
        }))
    return populated
