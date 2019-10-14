import os
import warnings


from Bio.Seq import Seq
import pytest



@pytest.fixture
def seq():
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


def test_three_frame_translation(seq):
    from khtools.partition_reads import three_frame_translation

    test = [str(x) for x in three_frame_translation(seq)]
    true = ['RLLNTDINNIRKIAI*L*ILFC', 'ACLILTSIILGKSQYNCKSCSV',
            'LA*Y*HQ*Y*ENRNITVNPVL']
    assert test == true


def test_three_frame_translation_no_stops(seq):
    from khtools.partition_reads import three_frame_translation_no_stops

    test = {k: str(v) for k, v in
            three_frame_translation_no_stops(seq).items()}
    true = {2: 'ACLILTSIILGKSQYNCKSCSV'}
    assert test == true


def test_six_frame_translation_no_stops(seq):
    from khtools.partition_reads import six_frame_translation_no_stops

    test = {k: str(v) for k, v in
            six_frame_translation_no_stops(seq).items()}
    true = {2: 'ACLILTSIILGKSQYNCKSCSV',
            -2: 'TEQDLQLYCDFPNIIDVSIKQA',
            -3: 'QNRIYSYIAIFLILLMSVLSK'}
    assert test == true

