"""
test_sequence_encodings.py

Tests for re-encoding biological sequence data
"""
import pytest

from sencha.sequence_encodings import VALID_PEPTIDE_MOLECULES

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


@pytest.fixture(params=VALID_PEPTIDE_MOLECULES)
def reduced_alphabet(request):
    return request.param


def test_translations():
    from sencha.sequence_encodings import (
        DAYHOFF_MAPPING,
        HP_MAPPING,
        BOTVINNIK_MAPPING,
        AMINO_ACID_SINGLE_LETTERS,
        DNA_ALPHABET,
        AMINO_KETO_MAPPING,
        PURINE_PYRIMIDINE_MAPPING,
        WEAK_STRONG_MAPPING,
    )

    assert all(x in DAYHOFF_MAPPING for x in AMINO_ACID_SINGLE_LETTERS)
    assert all(x in HP_MAPPING for x in AMINO_ACID_SINGLE_LETTERS)
    assert all(x in BOTVINNIK_MAPPING for x in AMINO_ACID_SINGLE_LETTERS)

    assert all(x in AMINO_KETO_MAPPING for x in DNA_ALPHABET)
    assert all(x in PURINE_PYRIMIDINE_MAPPING for x in DNA_ALPHABET)
    assert all(x in WEAK_STRONG_MAPPING for x in DNA_ALPHABET)


# -------------------- Test nucleotide encodings ---------------------------- #
def test_amino_keto_ize(nucleotide_string):
    from sencha.sequence_encodings import amino_keto_ize

    test = amino_keto_ize(nucleotide_string)
    true = "KMKKMMM"
    assert test == true


def test_weak_strong_ize(nucleotide_string):
    from sencha.sequence_encodings import weak_strong_ize

    test = weak_strong_ize(nucleotide_string)
    true = "SWWWWSW"
    assert test == true


def test_purine_pyrimidize(nucleotide_string):
    from sencha.sequence_encodings import purine_pyrimidize

    test = purine_pyrimidize(nucleotide_string)
    true = "RRYYRYR"
    assert test == true


# -------------------- Test peptide encodings ---------------------------- #


def test_dayhoffize(peptide_string):
    from sencha.sequence_encodings import dayhoffize

    test = dayhoffize(peptide_string)
    true = "bbbdbfecdac"
    assert test == true


def test_dayhoff_v2_ize(peptide_string):
    from sencha.sequence_encodings import dayhoff_v2_ize

    test = dayhoff_v2_ize(peptide_string)
    true = "BbBdbfecdac"
    assert test == true


def test_hpize(peptide_string):
    from sencha.sequence_encodings import hpize

    test = hpize(peptide_string)
    true = "phpphhhpppp"
    assert test == true


def test_botvinnikize(peptide_string):
    from sencha.sequence_encodings import botvinnikize

    test = botvinnikize(peptide_string)
    true = "dadkacbfghf"
    assert test == true


def test_peptide_constants():
    from sencha.sequence_encodings import PEPTIDE_MAPPINGS, AMINO_ACID_SINGLE_LETTERS

    for key, mapping in PEPTIDE_MAPPINGS.items():
        try:
            assert all(x in mapping for x in AMINO_ACID_SINGLE_LETTERS)
        except AssertionError:
            not_found = [x for x in AMINO_ACID_SINGLE_LETTERS if x not in mapping]
            raise AssertionError(
                f"Amino acids not present in {key} mapping: " f"{not_found}"
            )


def test_encode_peptide(peptide_string, reduced_alphabet):
    from sencha.sequence_encodings import encode_peptide

    test = encode_peptide(peptide_string, reduced_alphabet)
    true = peptide_string
    if reduced_alphabet == "dayhoff":
        true = "bbbdbfecdac"
    if reduced_alphabet == "dayhoff6":
        true = "bbbdbfecdac"
    elif reduced_alphabet == "hydrophobic-polar" or reduced_alphabet == "hp":
        true = "phpphhhpppp"
    elif reduced_alphabet == "hydrophobic-polar2" or reduced_alphabet == "hp2":
        true = "phpphhhpppp"
    elif reduced_alphabet == "protein" or reduced_alphabet == "peptide":
        true = peptide_string
    elif reduced_alphabet == "protein20" or reduced_alphabet == "aa20":
        true = peptide_string
    elif reduced_alphabet == "botvinnik":
        true = "dadkacbfghf"
    elif reduced_alphabet == "botvinnik8":
        true = "dadkacbfghf"
    elif reduced_alphabet == "aa9":
        true = "aaafabbdgbd"
    elif reduced_alphabet == "gbmr4":
        true = "aaababbaaba"
    elif reduced_alphabet == "sdm12":
        true = "eaejafgcchc"
    elif reduced_alphabet == "hsdm17":
        true = "gagqajkcdmc"
    assert test == true
