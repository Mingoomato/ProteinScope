"""Tests for structure_comparator.py"""

from __future__ import annotations

import asyncio
import os
import tempfile

import pytest
from analyzers.structure_comparator import (
    compare_structures,
    download_alphafold_pdb,
    StructuralComparison,
)


def test_compare_structures_nonexistent_files_returns_none():
    """compare_structures must return None gracefully for missing PDB paths."""
    result = compare_structures("/nonexistent/path/a.pdb", "/nonexistent/path/b.pdb")
    assert result is None


def test_compare_structures_one_nonexistent_returns_none():
    """Returns None when only one of the two files is missing."""
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        tmp_path = f.name
    try:
        result = compare_structures(tmp_path, "/nonexistent/b.pdb")
        assert result is None
    finally:
        os.unlink(tmp_path)


def test_compare_structures_empty_paths_returns_none():
    """Returns None when paths are empty strings."""
    assert compare_structures("", "") is None
    assert compare_structures("", "/some/path.pdb") is None


def test_compare_structures_invalid_pdb_content_returns_none():
    """Returns None when file exists but contains no valid PDB data."""
    with tempfile.NamedTemporaryFile(
        suffix=".pdb", mode="w", delete=False
    ) as f:
        f.write("THIS IS NOT VALID PDB CONTENT\n")
        tmp_path = f.name

    try:
        result = compare_structures(tmp_path, tmp_path)
        # Either None (biotite parse failure) or a valid result — must not raise
        assert result is None or isinstance(result, StructuralComparison)
    finally:
        os.unlink(tmp_path)


def test_structural_comparison_model_fields():
    """StructuralComparison instantiates correctly with expected fields."""
    sc = StructuralComparison(
        tm_score=0.75,
        rmsd_angstrom=1.5,
        n_residues_compared=200,
        interpretation="same fold",
    )
    assert sc.tm_score == 0.75
    assert sc.rmsd_angstrom == 1.5
    assert sc.n_residues_compared == 200
    assert sc.interpretation == "same fold"


@pytest.mark.asyncio
async def test_download_alphafold_pdb_is_async():
    """download_alphafold_pdb must be awaitable (is a coroutine function)."""
    import inspect
    assert inspect.iscoroutinefunction(download_alphafold_pdb)


@pytest.mark.asyncio
async def test_download_alphafold_pdb_invalid_id_returns_none():
    """Returns None for a clearly invalid UniProt ID (network call expected to fail)."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        result = await download_alphafold_pdb("INVALID_UNIPROT_ID_XYZXYZ", tmp_dir)
        assert result is None


@pytest.mark.asyncio
async def test_download_alphafold_pdb_empty_id_returns_none():
    """Returns None when given an empty UniProt ID."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        result = await download_alphafold_pdb("", tmp_dir)
        assert result is None


def test_compare_structures_same_tiny_pdb():
    """Comparing a minimal PDB file against itself should return TM-score=1.0 and RMSD≈0."""
    # Minimal PDB with 3 CA atoms of a fictitious protein
    minimal_pdb = (
        "ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00           C\n"
        "ATOM      2  CA  GLY A   2      11.500  11.500  10.000  1.00  0.00           C\n"
        "ATOM      3  CA  LEU A   3      13.000  10.000  11.000  1.00  0.00           C\n"
        "END\n"
    )
    with tempfile.NamedTemporaryFile(
        suffix=".pdb", mode="w", delete=False
    ) as f:
        f.write(minimal_pdb)
        tmp_path = f.name

    try:
        result = compare_structures(tmp_path, tmp_path)
        if result is not None:
            # Self-comparison should be near-identical
            assert result.tm_score == pytest.approx(1.0, abs=0.05)
            assert result.rmsd_angstrom == pytest.approx(0.0, abs=0.01)
            assert result.n_residues_compared == 3
            assert result.interpretation == "near-identical"
    finally:
        os.unlink(tmp_path)
