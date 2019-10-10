"""
test_os_utils.py

Tests for operating system utilities
"""

import os

import pytest


# Fixtures are functions-turned-variables that can be used across multiple
# tests. conftest.py contains fixtures that can be used by any test file
@pytest.fixture
def folder():
    return "test-folder"

# -------------------- Test nucleotide encodings ---------------------------- #
def test_amino_keto_ize():
    from khtools.sequence_encodings import amino_keto_ize

    test = amino_keto_ize("GATTACA")
    true = 'KMKKMMM'
    assert test == true

def test_weak_strong_ize():
    from khtools.sequence_encodings import weak_strong_ize

    test = weak_strong_ize("GATTACA")
    true = 'SWWWWSW'
    assert test == true

def test_purine_pyrimidize():
    from khtools.sequence_encodings import purine_pyrimidize

    test = purine_pyrimidize("GATTACA")
    true = 'RRYYRYR'
    assert test == true

# -------------------- Test peptide encodings ---------------------------- #
def test_dayhoffize():
    from khtools.sequence_encodings import dayhoffize

    test = dayhoffize("SASHAFIERCE")
    true = 'bbbdbfecdac'
    assert test == true


def test_dayhoff_v2_ize():
    from khtools.sequence_encodings import dayhoff_v2_ize

    test = dayhoff_v2_ize("SASHAFIERCE")
    true = 'BbBdbfecdac'
    assert test == true


def test_hpize():
    from khtools.sequence_encodings import hpize

    test = hpize("SASHAFIERCE")
    true = 'phpphhhpppp'
    assert test == true


def test_botvinnikize():
    from khtools.sequence_encodings import botvinnikize

    test = botvinnikize("SASHAFIERCE")
    true = 'dadkacbfghf'
    assert test == true

