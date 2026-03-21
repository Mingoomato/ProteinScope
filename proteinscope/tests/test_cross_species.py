"""Tests for analyzers/cross_species.py."""

from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from analyzers.cross_species import compute_identity


def test_compute_identity_identical():
    assert compute_identity("ACDEF", "ACDEF") == 100.0


def test_compute_identity_completely_different():
    assert compute_identity("AAAAA", "BBBBB") == 0.0


def test_compute_identity_partial():
    result = compute_identity("ACDEF", "ACXYZ")
    # First 2 of 5 match => 40%
    assert result == 40.0


def test_compute_identity_different_lengths():
    # Shorter seq compared against longer: max_len used as denominator
    result = compute_identity("ACD", "ACDEF")
    # 3 matches out of max length 5 => 60%
    assert result == 60.0


def test_compute_identity_empty():
    assert compute_identity("", "ABC") == 0.0
    assert compute_identity("ABC", "") == 0.0


def test_compute_identity_single_residue():
    assert compute_identity("A", "A") == 100.0
    assert compute_identity("A", "B") == 0.0
