"""
test_sequence_encodings.py

Tests for re-encoding biological sequence data
"""
import pytest


# Fixtures are functions-turned-variables that can be used across multiple
# tests. conftest.py contains fixtures that can be used by any test file
@pytest.fixture
def folder():
    return "test-folder"


@pytest.fixture
def peptide_string():
    return "SASHAFIERCE"


@pytest.fixture
def nucleotide_string():
    return "GATTACA"


def test_translations():
    from khtools.sequence_encodings import DAYHOFF_MAPPING, HP_MAPPING, \
        BOTVINNIK_MAPPING, AMINO_ACID_SINGLE_LETTERS, DNA_ALPHABET, \
        AMINO_KETO_MAPPING, PURINE_PYRIMIDINE_MAPPING, WEAK_STRONG_MAPPING

    assert all(x in DAYHOFF_MAPPING for x in AMINO_ACID_SINGLE_LETTERS)
    assert all(x in HP_MAPPING for x in AMINO_ACID_SINGLE_LETTERS)
    assert all(x in BOTVINNIK_MAPPING for x in AMINO_ACID_SINGLE_LETTERS)

    assert all(x in AMINO_KETO_MAPPING for x in DNA_ALPHABET)
    assert all(x in PURINE_PYRIMIDINE_MAPPING for x in DNA_ALPHABET)
    assert all(x in WEAK_STRONG_MAPPING for x in DNA_ALPHABET)


# -------------------- Test nucleotide encodings ---------------------------- #
def test_amino_keto_ize(nucleotide_string):
    from khtools.sequence_encodings import amino_keto_ize

    test = amino_keto_ize(nucleotide_string)
    true = 'KMKKMMM'
    assert test == true

def test_weak_strong_ize(nucleotide_string):
    from khtools.sequence_encodings import weak_strong_ize

    test = weak_strong_ize(nucleotide_string)
    true = 'SWWWWSW'
    assert test == true

def test_purine_pyrimidize(nucleotide_string):
    from khtools.sequence_encodings import purine_pyrimidize

    test = purine_pyrimidize(nucleotide_string)
    true = 'RRYYRYR'
    assert test == true

# -------------------- Test peptide encodings ---------------------------- #
def test_dayhoffize(peptide_string):
    from khtools.sequence_encodings import dayhoffize

    test = dayhoffize(peptide_string)
    true = 'bbbdbfecdac'
    assert test == true


def test_dayhoff_v2_ize(peptide_string):
    from khtools.sequence_encodings import dayhoff_v2_ize

    test = dayhoff_v2_ize(peptide_string)
    true = 'BbBdbfecdac'
    assert test == true


def test_hpize(peptide_string):
    from khtools.sequence_encodings import hpize

    test = hpize(peptide_string)
    true = 'phpphhhpppp'
    assert test == true


def test_botvinnikize(peptide_string):
    from khtools.sequence_encodings import botvinnikize

    test = botvinnikize(peptide_string)
    true = 'dadkacbfghf'
    assert test == true


def test_encode_peptide(peptide_string, molecule):
    from khtools.sequence_encodings import encode_peptide

    test = encode_peptide(peptide_string, molecule)
    if molecule == 'dayhoff':
        true = 'bbbdbfecdac'
    elif molecule == 'hydrophobic-polar':
        true = 'phpphhhpppp'
    elif molecule == 'protein':
        true = peptide_string
    assert test == true

