"""Domain architecture scanner using pyhmmer + Pfam-A HMM profiles.

pyhmmer 0.12.0 is installed. Pfam-A.hmm must be downloaded separately to
proteinscope/.data/pfam/Pfam-A.hmm (~500 MB).

If the database is not present, falls back to returning an empty list with
a descriptive message. Call ``download_pfam_info()`` to get fetch instructions.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field

# ---------------------------------------------------------------------------
# Pfam database location — relative to this file's package root
# ---------------------------------------------------------------------------
_PFAM_RELATIVE = Path(__file__).resolve().parent.parent / ".data" / "pfam" / "Pfam-A.hmm"


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class DomainHit(BaseModel):
    name: str
    accession: str
    start: int
    end: int
    evalue: float
    score: float
    coverage: float
    description: str = ""


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def get_pfam_path() -> Optional[Path]:
    """Return the absolute path to Pfam-A.hmm, or None if it does not exist.

    Checks the canonical location at ``<project_root>/.data/pfam/Pfam-A.hmm``.
    Set the environment variable ``PFAM_HMM_PATH`` to override the default.
    """
    env_override = os.environ.get("PFAM_HMM_PATH")
    if env_override:
        p = Path(env_override)
        if p.is_file():
            return p

    if _PFAM_RELATIVE.is_file():
        return _PFAM_RELATIVE

    return None


def download_pfam_info() -> str:
    """Return human-readable instructions for downloading Pfam-A.hmm.

    The Pfam database is ~500 MB and must be downloaded once before domain
    scanning is available.
    """
    target = str(_PFAM_RELATIVE)
    return (
        "Pfam-A.hmm not found. To enable domain scanning, download it:\n\n"
        "  # Option 1: wget / curl (Linux/macOS)\n"
        "  wget -O Pfam-A.hmm.gz "
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz\n"
        "  gunzip Pfam-A.hmm.gz\n\n"
        "  # Option 2: Python\n"
        "  import urllib.request, gzip, shutil\n"
        "  url = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'\n"
        "  urllib.request.urlretrieve(url, 'Pfam-A.hmm.gz')\n"
        "  with gzip.open('Pfam-A.hmm.gz', 'rb') as f_in, open('Pfam-A.hmm', 'wb') as f_out:\n"
        "      shutil.copyfileobj(f_in, f_out)\n\n"
        f"  # Then move it to: {target}\n\n"
        "Alternatively, set the PFAM_HMM_PATH environment variable to the "
        "absolute path of your existing Pfam-A.hmm file."
    )


# ---------------------------------------------------------------------------
# Core scanning function
# ---------------------------------------------------------------------------

def scan_pfam_domains(sequence: str) -> list[DomainHit]:
    """Scan a protein sequence against Pfam-A using pyhmmer.

    Requires Pfam-A.hmm at the canonical location (see ``get_pfam_path()``).
    Returns an empty list if the database is absent or scanning fails.

    Args:
        sequence: Single-letter amino acid sequence string.

    Returns:
        List of DomainHit objects sorted by E-value (best first),
        filtered to E-value < 0.001.
    """
    if not sequence:
        return []

    pfam_path = get_pfam_path()
    if pfam_path is None:
        return []

    try:
        import pyhmmer
        from pyhmmer.easel import TextSequence, Alphabet
        from pyhmmer.plan7 import HMMFile

        alphabet = Alphabet.amino()
        text_seq = TextSequence(name=b"query", sequence=sequence.encode())
        digital_seq = text_seq.digitize(alphabet)

        hits_out: list[DomainHit] = []

        with HMMFile(str(pfam_path)) as hmm_file:
            for hits in pyhmmer.hmmscan(digital_seq, hmm_file, cpus=1):
                for hit in hits:
                    if hit.evalue >= 0.001:
                        continue

                    # Extract domain coordinates from the best domain
                    for domain in hit.domains:
                        if not domain.included:
                            continue

                        start = int(domain.env_from)
                        end = int(domain.env_to)
                        domain_len = end - start + 1
                        coverage = domain_len / len(sequence) if len(sequence) > 0 else 0.0

                        accession = hit.accession.decode() if hit.accession else ""
                        name = hit.name.decode() if hit.name else ""
                        description = hit.description.decode() if hit.description else ""

                        hits_out.append(DomainHit(
                            name=name,
                            accession=accession,
                            start=start,
                            end=end,
                            evalue=float(hit.evalue),
                            score=float(hit.score),
                            coverage=round(coverage, 4),
                            description=description,
                        ))

        hits_out.sort(key=lambda h: h.evalue)
        return hits_out

    except ImportError:
        # pyhmmer not installed
        return []
    except Exception:
        return []


# ---------------------------------------------------------------------------
# Architecture comparison
# ---------------------------------------------------------------------------

def compare_domain_architectures(
    domains_a: list[DomainHit],
    domains_b: list[DomainHit],
) -> dict:
    """Compare the domain architectures of two proteins.

    Matches domains by Pfam accession (e.g. PF00001).

    Args:
        domains_a: Domain hits for protein A.
        domains_b: Domain hits for protein B.

    Returns:
        Dict with keys ``shared``, ``only_in_a``, ``only_in_b``, each a list
        of accession strings.
    """
    acc_a: set[str] = {d.accession for d in domains_a if d.accession}
    acc_b: set[str] = {d.accession for d in domains_b if d.accession}

    shared = sorted(acc_a & acc_b)
    only_in_a = sorted(acc_a - acc_b)
    only_in_b = sorted(acc_b - acc_a)

    return {
        "shared": shared,
        "only_in_a": only_in_a,
        "only_in_b": only_in_b,
    }
