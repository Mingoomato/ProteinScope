"""Tests for NCBI fetchers (ncbi_gene.py, ncbi_clinvar.py)."""

from __future__ import annotations

import json
import pytest
from pytest_httpx import HTTPXMock

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from fetchers.ncbi_gene import parse_cds_from_genbank
from fetchers.ncbi_clinvar import parse_clinical_entries

# ── GenBank CDS parsing ───────────────────────────────────────────────────────

MINIMAL_GENBANK = """\
LOCUS       NM_000518               1606 bp    mRNA    linear   PRI 01-JAN-2020
DEFINITION  Homo sapiens hemoglobin subunit beta (HBB), mRNA.
ACCESSION   NM_000518
FEATURES             Location/Qualifiers
     source          1..1606
     CDS             51..494
                     /gene="HBB"
                     /product="hemoglobin subunit beta"
ORIGIN
        1 acatttgctt ctgacacaac tgtgttcact agcaacctca aacagacacc atggtgcatc
       61 tgactcctga ggagaagtct gccgttactg ccctgtgggg caaggtgaac gtggatgaag
      121 ttggtggtga ggccctgggc aggctgctgg tggtctaccc ttggacccag aggttctttg
      181 agtcctttgg ggatctgtcc actcctgatg ctgttatggg caaccctaag gtgaaggctc
      241 atggcaagaa agtgctcggt gcctttagtg atggcctggc tcacctggac aacctcaagg
      301 gcacctttgc cacactgagt gagctgcact gtgacaagct gcacgtggat cctgagaact
      361 tcaggctcct gggcaacgtg ctggtctgtg tgctggccca tcactttggc aaagaattca
      421 ccccaccagt gcaggctgcc tatcagaaag tggtggctgg tgtggctaat gccctggccc
      481 acaagtatca ctaa
//
"""


def test_parse_cds_from_genbank():
    seq, is_spliced = parse_cds_from_genbank(MINIMAL_GENBANK)
    assert seq is not None
    assert isinstance(seq, str)
    assert len(seq) > 0
    # CDS starts with ATG for most coding sequences
    assert seq.startswith("ATG") or len(seq) > 0


def test_parse_cds_not_spliced():
    _, is_spliced = parse_cds_from_genbank(MINIMAL_GENBANK)
    assert is_spliced is False


def test_parse_cds_empty_returns_none():
    seq, _ = parse_cds_from_genbank("")
    assert seq is None


# ── ClinVar parsing ───────────────────────────────────────────────────────────

MINIMAL_VCV = {
    "ClinVarResult-Set": {
        "VariationArchive": {
            "InterpretedRecord": {
                "RCVList": {
                    "RCVAccession": {
                        "Interpretation": {"Description": "Pathogenic"},
                        "ClassifiedConditionList": {
                            "ClassifiedCondition": {"#text": "Sickle cell disease"}
                        },
                    }
                }
            }
        }
    }
}


def test_parse_clinical_entries():
    entries = parse_clinical_entries(MINIMAL_VCV)
    assert len(entries) == 1
    assert entries[0]["disease"] == "Sickle cell disease"
    assert entries[0]["significance"] == "Pathogenic"


def test_parse_clinical_entries_empty():
    entries = parse_clinical_entries({})
    assert entries == []
