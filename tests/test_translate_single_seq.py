import pytest
from sencha.translate_single_seq import TranslateSingleSeq


class TestTranslateSingleSeq:
    def __init__(self, seq, verbose):
        self.translate_ss = TranslateSingleSeq(seq, verbose)


@pytest.fixture()
def translation(seq):
    return TestTranslateSingleSeq(seq, True)


def test__reverse_complement(translation):
    # Palindromic sequence
    assert translation.translate_ss._reverse_complement("AACCGGTT") == "AACCGGTT"

    # Non-palindrome
    assert translation.translate_ss._reverse_complement("AACC") == "GGTT"


def test__translate_codon(translation):
    # Make sure first base as "N" results in "X ambiguous amino acid
    assert translation.translate_ss._translate_codon("NAA") == "X"

    # Threonine can have ambiguous final base
    assert translation.translate_ss._translate_codon("ACN") == "T"

    # Leucine can have an ambiguous final base
    assert translation.translate_ss._translate_codon("CTN") == "L"

    # Valine can have ambiguous final base
    assert translation.translate_ss._translate_codon("GTN") == "V"

    # Alanine can have ambiguous final base
    assert translation.translate_ss._translate_codon("GCN") == "A"

    # Arginine can have ambiguous final base
    assert translation.translate_ss._translate_codon("CGN") == "R"

    # Glycine can have ambiguous final base
    assert translation.translate_ss._translate_codon("GGN") == "G"

    # Make sure too short codons don't translate
    assert translation.translate_ss._translate_codon("AC") == ""


def test_three_frame_translation(translation):
    test = [str(x) for x in translation.translate_ss.three_frame_translation()]
    true = ["RLLNTDINNIRKIAI*L*ILFC", "ACLILTSIILGKSQYNCKSCSV", "LA*Y*HQ*Y*ENRNITVNPVL"]
    assert test == true


def test_three_frame_fwd_translation(translation):
    test = {
        k: str(v)
        for k, v in translation.translate_ss.three_frame_translation_no_stops(1).items()
    }
    true = {2: "ACLILTSIILGKSQYNCKSCSV"}
    assert test == true


def test_three_frame_rev_translation(translation):
    result = translation.translate_ss.three_frame_translation_no_stops(-1)
    test = {k: str(v) for k, v in result.items()}
    true = {-2: "TEQDLQLYCDFPNIIDVSIKQA", -3: "QNRIYSYIAIFLILLMSVLSK"}
    assert test == true


def test_six_frame_translation_no_stops(translation):
    result = translation.translate_ss.six_frame_translation_no_stops()
    test = {k: str(v) for k, v in result.items()}
    true = {
        2: "ACLILTSIILGKSQYNCKSCSV",
        -2: "TEQDLQLYCDFPNIIDVSIKQA",
        -3: "QNRIYSYIAIFLILLMSVLSK",
    }
    assert test == true


def test_three_frame_translation_stops(translation):
    result = translation.translate_ss.three_frame_translation_stops(-1)
    test = {k: str(v) for k, v in result.items()}
    true = {
        -1: "DRTGFTVILRFS*YY*CQY*AS",
        -2: "TEQDLQLYCDFPNIIDVSIKQA",
        -3: "QNRIYSYIAIFLILLMSVLSK",
    }
    assert test == true


def test_six_frame_translation_stops(translation):
    result = translation.translate_ss.six_frame_translation()
    test = {k: str(v) for k, v in result.items()}
    true = {
        1: "RLLNTDINNIRKIAI*L*ILFC",
        2: "ACLILTSIILGKSQYNCKSCSV",
        3: "LA*Y*HQ*Y*ENRNITVNPVL",
        -1: "DRTGFTVILRFS*YY*CQY*AS",
        -2: "TEQDLQLYCDFPNIIDVSIKQA",
        -3: "QNRIYSYIAIFLILLMSVLSK",
    }
    assert test == true
