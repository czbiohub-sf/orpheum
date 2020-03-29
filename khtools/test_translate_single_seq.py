import warnings
import pytest
from khtools.translate_single_seq import TranslateSingleSeq


@pytest.fixture
def seq():
    from Bio.Seq import Seq
    s = 'CGCTTGCTTAATACTGACATCAATAATATTAGGAAAATCGCAATATAACTGTAAATCCTGTTCTGTC'
    with warnings.catch_warnings():
        # Ignore The following warning because we don't use Bio.Alphabet
        # explicitly:
        # PendingDeprecationWarning: We intend to remove or replace
        # Bio.Alphabet in 2020, ideally avoid using it explicitly in your
        # code. Please get in touch if you will be adversely affected by this.
        # https://github.com/biopython/biopython/issues/2046
        warnings.simplefilter("ignore")
        return Seq(s)


class TestTranslateSingleSeq:
    def __init__(self, seq, verbose):
        self.translate_ss = TranslateSingleSeq(seq, verbose)


@pytest.fixture()
def translation(seq):
    return TestTranslateSingleSeq(seq, True)


def test_three_frame_translation(translation):
    test = [str(x) for x in translation.translate_ss.three_frame_translation()]
    true = [
        'RLLNTDINNIRKIAI*L*ILFC', 'ACLILTSIILGKSQYNCKSCSV',
        'LA*Y*HQ*Y*ENRNITVNPVL'
    ]
    assert test == true


def test_three_frame_fwd_translation(translation):
    test = {
        k: str(v)
        for k, v in translation.translate_ss.three_frame_translation_no_stops(
            1).items()
    }
    true = {2: 'ACLILTSIILGKSQYNCKSCSV'}
    assert test == true


def test_three_frame_rev_translation(translation):
    result = translation.translate_ss.three_frame_translation_no_stops(-1)
    test = {
        k: str(
            v) for k, v in result.items()}
    true = {
        -2: 'TEQDLQLYCDFPNIIDVSIKQA',
        -3: 'QNRIYSYIAIFLILLMSVLSK'
    }
    assert test == true


def test_six_frame_translation(translation):
    result = translation.translate_ss.six_frame_translation_no_stops()
    test = {
        k: str(
            v) for k, v in result.items()}
    true = {
        2: 'ACLILTSIILGKSQYNCKSCSV',
        -2: 'TEQDLQLYCDFPNIIDVSIKQA',
        -3: 'QNRIYSYIAIFLILLMSVLSK'
    }
    assert test == true
