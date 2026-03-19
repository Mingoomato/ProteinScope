"""Sequence alignment engine — pairwise and multiple sequence alignment.

Uses Bio.Align.PairwiseAligner (BioPython 1.80+) for pairwise alignment.
MSA uses MAFFT subprocess when available, falls back to BioPython pairwise.
All functions return None/empty gracefully if alignment fails.
"""

from __future__ import annotations

import math
import subprocess
from typing import Optional

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Physicochemical conservation groups for similarity scoring
# ---------------------------------------------------------------------------
_CONSERVATIVE_GROUPS: list[frozenset[str]] = [
    frozenset("LIFVM"),
    frozenset("HKR"),
    frozenset("ED"),
    frozenset("NQST"),
    frozenset("AG"),
    frozenset("CMPW"),
    frozenset("FYW"),
]


def _is_conservative(a: str, b: str) -> bool:
    """Return True if residues a and b fall in the same physicochemical group."""
    for group in _CONSERVATIVE_GROUPS:
        if a in group and b in group:
            return True
    return False


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class SiteConservation(BaseModel):
    site_description: str
    residue_a: str
    residue_b: str
    conserved: bool
    conservative: bool


class AlignmentResult(BaseModel):
    identity_pct: float
    similarity_pct: float
    gap_pct: float
    aligned_query: str
    aligned_target: str
    score: float
    site_conservation: list[SiteConservation] = Field(default_factory=list)


class MSAResult(BaseModel):
    aligned_sequences: dict[str, str]
    conservation_scores: list[float]  # per column, 1.0 = fully conserved
    n_sequences: int


# ---------------------------------------------------------------------------
# PairwiseAligner — initialised once at module level
# ---------------------------------------------------------------------------

def _make_aligner():
    """Build a BLOSUM62 global PairwiseAligner, returning None on import error."""
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.mode = "global"
        return aligner
    except Exception:
        return None


_ALIGNER = _make_aligner()


def _make_nucleotide_aligner():
    """Build a global PairwiseAligner for DNA/RNA using NUC44 (or match/mismatch fallback)."""
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
        aligner = PairwiseAligner()
        try:
            aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")
        except Exception:
            # NUC44 not available — use simple match/mismatch scoring
            aligner.match_score = 2.0
            aligner.mismatch_score = -1.0
        aligner.open_gap_score = -2.0
        aligner.extend_gap_score = -0.5
        aligner.mode = "global"
        return aligner
    except Exception:
        return None


_NUCLEOTIDE_ALIGNER = _make_nucleotide_aligner()

# Nucleotide characters (IUPAC including degenerate bases + RNA U)
_NUCL_CHARS = frozenset("ACGTUacgtuRYSWKMBDHVNryswkmbdhvn")


def _is_nucleotide_sequence(seq: str) -> bool:
    """Return True if >90% of non-gap characters are IUPAC nucleotide codes."""
    chars = [c for c in seq.upper() if c not in ("-", " ", "\n")]
    if len(chars) < 5:
        return False
    nucl_count = sum(1 for c in chars if c in _NUCL_CHARS)
    return nucl_count / len(chars) >= 0.90


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def align_pairwise(seq_a: str, seq_b: str) -> Optional[AlignmentResult]:
    """Perform global pairwise alignment between two protein sequences.

    Uses Bio.Align.PairwiseAligner with BLOSUM62 substitution matrix.

    Args:
        seq_a: First protein sequence (single-letter amino acid codes).
        seq_b: Second protein sequence.

    Returns:
        AlignmentResult with identity/similarity/gap percentages and aligned
        strings, or None if either sequence is empty or alignment fails.
    """
    if not seq_a or not seq_b:
        return None

    if _ALIGNER is None:
        return None

    try:
        alignments = _ALIGNER.align(seq_a, seq_b)
        # Take the best (first) alignment
        best = next(iter(alignments), None)
        if best is None:
            return None

        aligned_q = str(best[0])
        aligned_t = str(best[1])

        align_len = len(aligned_q)
        if align_len == 0:
            return None

        matches = 0
        similar = 0
        gaps = 0

        for rq, rt in zip(aligned_q, aligned_t):
            if rq == "-" or rt == "-":
                gaps += 1
            elif rq == rt:
                matches += 1
                similar += 1
            elif _is_conservative(rq, rt):
                similar += 1

        identity_pct = 100.0 * matches / align_len
        similarity_pct = 100.0 * similar / align_len
        gap_pct = 100.0 * gaps / align_len

        return AlignmentResult(
            identity_pct=round(identity_pct, 2),
            similarity_pct=round(similarity_pct, 2),
            gap_pct=round(gap_pct, 2),
            aligned_query=aligned_q,
            aligned_target=aligned_t,
            score=float(best.score),
        )

    except Exception:
        return None


def align_pairwise_nucleotide(seq_a: str, seq_b: str) -> Optional[AlignmentResult]:
    """Perform global pairwise alignment between two DNA or RNA sequences.

    Uses NUC44 matrix when available, otherwise match(+2)/mismatch(-1) scoring.
    T and U are treated as equivalent (both map to T for comparison).

    Args:
        seq_a: First nucleotide sequence (DNA or RNA, IUPAC codes).
        seq_b: Second nucleotide sequence.

    Returns:
        AlignmentResult with identity/gap percentages and gap-aligned strings,
        or None on failure.
    """
    if not seq_a or not seq_b:
        return None
    if _NUCLEOTIDE_ALIGNER is None:
        return None

    # Normalise RNA → DNA (U → T) for alignment, then restore in output
    def _norm(s: str) -> str:
        return s.upper().replace("U", "T")

    norm_a = _norm(seq_a)
    norm_b = _norm(seq_b)

    try:
        alignments = _NUCLEOTIDE_ALIGNER.align(norm_a, norm_b)
        best = next(iter(alignments), None)
        if best is None:
            return None

        aligned_q = str(best[0])
        aligned_t = str(best[1])

        # Restore original case/U if the input was RNA
        def _restore(aligned: str, original: str) -> str:
            is_rna = "U" in original.upper()
            result, orig_idx = [], 0
            for ch in aligned:
                if ch == "-":
                    result.append("-")
                else:
                    orig_ch = original[orig_idx] if orig_idx < len(original) else ch
                    result.append(orig_ch)
                    orig_idx += 1
            return "".join(result)

        aligned_q = _restore(aligned_q, seq_a)
        aligned_t = _restore(aligned_t, seq_b)

        align_len = len(aligned_q)
        if align_len == 0:
            return None

        matches = gaps = 0
        for qa, ta in zip(aligned_q, aligned_t):
            if qa == "-" or ta == "-":
                gaps += 1
            elif qa.upper().replace("U", "T") == ta.upper().replace("U", "T"):
                matches += 1

        identity_pct = 100.0 * matches / align_len
        gap_pct = 100.0 * gaps / align_len

        return AlignmentResult(
            identity_pct=round(identity_pct, 2),
            similarity_pct=round(identity_pct, 2),  # for nucleotides, identity == similarity
            gap_pct=round(gap_pct, 2),
            aligned_query=aligned_q,
            aligned_target=aligned_t,
            score=float(best.score),
        )

    except Exception:
        return None


def align_pairwise_auto(seq_a: str, seq_b: str) -> Optional[AlignmentResult]:
    """Auto-detect sequence type and delegate to the correct aligner.

    If both sequences look like nucleotides (>90% IUPAC nt codes) uses
    align_pairwise_nucleotide; otherwise uses align_pairwise (protein/BLOSUM62).
    """
    if _is_nucleotide_sequence(seq_a) and _is_nucleotide_sequence(seq_b):
        return align_pairwise_nucleotide(seq_a, seq_b)
    return align_pairwise(seq_a, seq_b)


def align_at_annotated_sites(
    seq_a: str,
    seq_b: str,
    sites_a: list[dict],
) -> list[dict]:
    """Compare annotated functional sites between two aligned sequences.

    For each site in sites_a, finds the corresponding residues in seq_b via
    pairwise alignment and reports conservation at that site.

    Args:
        seq_a: Query sequence.
        seq_b: Target sequence to compare against.
        sites_a: List of site dicts with keys: start (int), end (int),
                 description (str). Positions are 0-based, end exclusive.

    Returns:
        List of dicts with site conservation information. Empty list on failure.
    """
    if not seq_a or not seq_b or not sites_a:
        return []

    result = align_pairwise(seq_a, seq_b)
    if result is None:
        return []

    aligned_q = result.aligned_query
    aligned_t = result.aligned_target

    # Build mapping from original seq_a positions → aligned position
    orig_to_aln: dict[int, int] = {}
    orig_idx = 0
    for aln_idx, ch in enumerate(aligned_q):
        if ch != "-":
            orig_to_aln[orig_idx] = aln_idx
            orig_idx += 1

    output: list[dict] = []

    for site in sites_a:
        try:
            start = int(site.get("start", 0))
            end = int(site.get("end", start + 1))
            description = str(site.get("description", ""))

            # For each position in the site window
            for pos_a in range(start, end):
                aln_pos = orig_to_aln.get(pos_a)
                if aln_pos is None:
                    continue

                res_a = aligned_q[aln_pos] if aln_pos < len(aligned_q) else "-"
                res_b = aligned_t[aln_pos] if aln_pos < len(aligned_t) else "-"

                conserved = res_a == res_b and res_a != "-"
                conservative = _is_conservative(res_a, res_b) if (res_a != "-" and res_b != "-") else False

                output.append({
                    "site_description": description,
                    "pos_a": pos_a,
                    "pos_b": _aln_to_orig_b(aligned_t, aln_pos),
                    "residue_a": res_a,
                    "residue_b": res_b,
                    "conserved": conserved,
                    "conservative": conservative,
                })

        except Exception:
            continue

    return output


def _aln_to_orig_b(aligned_t: str, aln_pos: int) -> int:
    """Convert aligned position to original index in sequence B."""
    count = 0
    for i, ch in enumerate(aligned_t):
        if i == aln_pos:
            return count if ch != "-" else -1
        if ch != "-":
            count += 1
    return -1


def build_msa(sequences: dict[str, str]) -> Optional[MSAResult]:
    """Build a multiple sequence alignment from a dict of named sequences.

    Tries MAFFT subprocess first (``mafft --quiet --auto --reorder``).
    Falls back to progressive pairwise alignment against the first sequence
    when MAFFT is unavailable (common on Windows).

    Args:
        sequences: Mapping of sequence name → sequence string.

    Returns:
        MSAResult with aligned sequences and per-column conservation scores
        (1.0 = fully conserved, 0.0 = maximally variable), or None on failure.
    """
    if not sequences or len(sequences) < 2:
        return None

    aligned: Optional[dict[str, str]] = _msa_mafft(sequences)
    if aligned is None:
        aligned = _msa_pairwise_fallback(sequences)
    if aligned is None:
        return None

    conservation = _compute_conservation(aligned)

    return MSAResult(
        aligned_sequences=aligned,
        conservation_scores=conservation,
        n_sequences=len(aligned),
    )


# ---------------------------------------------------------------------------
# Internal MSA helpers
# ---------------------------------------------------------------------------

def _msa_mafft(sequences: dict[str, str]) -> Optional[dict[str, str]]:
    """Run MAFFT on the input sequences. Returns None if MAFFT unavailable."""
    try:
        fasta_input = "".join(
            f">{name}\n{seq}\n" for name, seq in sequences.items()
        )
        proc = subprocess.run(
            ["mafft", "--quiet", "--auto", "--reorder", "-"],
            input=fasta_input,
            capture_output=True,
            text=True,
            timeout=60,
        )
        if proc.returncode != 0:
            return None

        return _parse_fasta(proc.stdout)

    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return None
    except Exception:
        return None


def _msa_pairwise_fallback(sequences: dict[str, str]) -> Optional[dict[str, str]]:
    """Pairwise progressive alignment: align all sequences to the first one.

    This is a simple approximation of MSA suitable as a fallback when MAFFT
    is not available. The first sequence acts as the reference/anchor.
    """
    names = list(sequences.keys())
    seqs = list(sequences.values())

    if not seqs:
        return None

    ref_name = names[0]
    ref_seq = seqs[0]

    aligned: dict[str, str] = {}

    # Align all sequences against the first and collect aligned strings
    ref_versions: list[str] = []
    other_aligned: list[tuple[str, str]] = []

    for name, seq in zip(names[1:], seqs[1:]):
        res = align_pairwise(ref_seq, seq)
        if res is None:
            # Pad shorter / longer sequences naively
            max_len = max(len(ref_seq), len(seq))
            aligned[name] = seq.ljust(max_len, "-")[:max_len]
            ref_versions.append(ref_seq.ljust(max_len, "-")[:max_len])
        else:
            ref_versions.append(res.aligned_query)
            other_aligned.append((name, res.aligned_target))

    if not ref_versions:
        aligned[ref_name] = ref_seq
        return aligned

    # Determine unified alignment length by merging gap positions
    # Build a consensus aligned reference by taking the longest ref alignment
    longest_ref = max(ref_versions, key=len)
    aln_len = len(longest_ref)

    aligned[ref_name] = longest_ref

    for (name, aln_seq), ref_v in zip(other_aligned, ref_versions):
        # Extend or trim to unified length
        if len(aln_seq) < aln_len:
            aln_seq = aln_seq + "-" * (aln_len - len(aln_seq))
        else:
            aln_seq = aln_seq[:aln_len]
        aligned[name] = aln_seq

    return aligned


def _parse_fasta(fasta_text: str) -> dict[str, str]:
    """Parse FASTA text into a dict of name → sequence."""
    result: dict[str, str] = {}
    current_name: Optional[str] = None
    buf: list[str] = []

    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None:
                result[current_name] = "".join(buf)
            current_name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)

    if current_name is not None:
        result[current_name] = "".join(buf)

    return result


def _compute_conservation(aligned: dict[str, str]) -> list[float]:
    """Compute per-column conservation as inverted Shannon entropy.

    Returns a list of floats in [0.0, 1.0] where 1.0 = fully conserved.
    """
    if not aligned:
        return []

    seqs = list(aligned.values())
    if not seqs:
        return []

    aln_len = max(len(s) for s in seqs)
    n_seqs = len(seqs)
    scores: list[float] = []

    for col in range(aln_len):
        column = [s[col] if col < len(s) else "-" for s in seqs]
        # Count residue frequencies (ignore gaps for entropy)
        residues = [r for r in column if r != "-"]
        if not residues:
            scores.append(0.0)
            continue

        counts: dict[str, int] = {}
        for r in residues:
            counts[r] = counts.get(r, 0) + 1

        n = len(residues)
        entropy = 0.0
        for cnt in counts.values():
            p = cnt / n
            if p > 0:
                entropy -= p * math.log2(p)

        max_entropy = math.log2(min(n, 20)) if n > 1 else 0.0
        if max_entropy == 0.0:
            conservation = 1.0
        else:
            conservation = max(0.0, 1.0 - entropy / max_entropy)

        scores.append(round(conservation, 4))

    return scores
