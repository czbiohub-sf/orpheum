from io import StringIO
import os
import warnings

from Bio.Seq import Seq
from click.testing import CliRunner
from khmer import Nodegraph
import pandas as pd
import pandas.util.testing as pdt
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


@pytest.fixture
def peptide_graph(data_folder):
    filename = os.path.join(data_folder, 'bloom_filter',
                            'Homo_sapiens.GRCh38.pep.subset.molecule-protein_ksize-7.bloomfilter.nodegraph')
    return Nodegraph.load(filename)


@pytest.fixture
def reads(data_folder):
    return os.path.join(data_folder,
                        'SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_first20.fq.gz')


def test_score_reads(reads, peptide_graph):
    from khtools.partition_reads import score_reads

    test = score_reads(reads, peptide_graph, peptide_ksize=7,
                                jaccard_threshold=0.9, molecule='protein')
    s = '''read_id,jaccard_in_peptide_db,n_kmers,classification
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1,1.0,16,coding
SRR306838.6196593 Ibis_Run100924_C3PO:6:29:16733:12435/1,0.0,17,non-coding
SRR306838.20767303 Ibis_Run100924_C3PO:6:104:6864:5062/1,0.0,16,non-coding
SRR306838.12582274 Ibis_Run100924_C3PO:6:62:11779:17975/1,0.0,16,non-coding
SRR306838.13334230 Ibis_Run100924_C3PO:6:66:16579:20350/1,0.0,12,non-coding
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1,1.0,18,coding
SRR306838.6813354 Ibis_Run100924_C3PO:6:32:10591:13073/1,0.0,17,non-coding
SRR306838.23113368 Ibis_Run100924_C3PO:6:114:13840:18459/1,0.0,14,non-coding
SRR306838.10872941 Ibis_Run100924_C3PO:6:53:6164:10522/1,0.0,16,non-coding
SRR306838.6192120 Ibis_Run100924_C3PO:6:29:5833:11991/1,0.0,16,non-coding
SRR306838.21295280 Ibis_Run100924_C3PO:6:106:2590:13965/1,0.0,17,non-coding
SRR306838.21201208 Ibis_Run100924_C3PO:6:106:2763:5109/1,0.0,0,non-coding
SRR306838.18327923 Ibis_Run100924_C3PO:6:92:9077:13885/1,0.0,16,non-coding
SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1,1.0,16,coding
SRR306838.21417895 Ibis_Run100924_C3PO:6:107:8793:5012/1,0.0,16,non-coding
SRR306838.17165743 Ibis_Run100924_C3PO:6:86:18789:18450/1,0.0,0,non-coding
SRR306838.21229494 Ibis_Run100924_C3PO:6:106:6163:7753/1,0.0,0,non-coding
SRR306838.21218773 Ibis_Run100924_C3PO:6:106:16921:6743/1,0.0,5,non-coding
SRR306838.20124664 Ibis_Run100924_C3PO:6:101:4701:5309/1,0.05555555555555555,18,non-coding
SRR306838.16841308 Ibis_Run100924_C3PO:6:85:6205:5805/1,0.0,17,non-coding'''
    true = pd.read_csv(StringIO(s))
    pdt.assert_equal(test, true)


def test_cli():

    runner = CliRunner()
    result = runner.invoke(hello, ['Peter'])
    assert result.exit_code == 0
    assert result.output == 'Hello Peter!\n'
