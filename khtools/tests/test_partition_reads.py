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
def reads(data_folder):
    return os.path.join(data_folder,
                        'SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_n22.fq.gz')


@pytest.fixture
def true_scores(data_folder, molecule, peptide_ksize):
    filename = os.path.join(data_folder, "partition_reads",
                            "SRR306838_GSM752691_hsa_br_F_1_trimmed_"
                            f"subsampled_n22__molecule-{molecule}_ksize-"
                            f"{peptide_ksize}.csv")
    return pd.read_csv(filename, index_col=0)


@pytest.fixture
def true_protein_coding_fasta():
    return """>SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2
TEQDLQLYCDFPNIIDVSIKQA
>SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1
QSSSPEFRVQSFSERTNARKKNNH
>SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2
LDPPYSRVITQRETENNQMTSE
"""


def test_score_reads(capsys, reads, peptide_graph, molecule, peptide_ksize,
                     true_scores, true_protein_coding_fasta):
    from khtools.partition_reads import score_reads

    test = score_reads(reads, peptide_graph, peptide_ksize=peptide_ksize,
                       molecule=molecule)
    pdt.assert_equal(test, true_scores)
    captured = capsys.readouterr()

    # Check that the proper sequences were output
    assert captured.out == true_protein_coding_fasta

    # Check tqdm iterations
    assert '22it' in captured.err


def test_cli(reads, peptide_bloom_filter, molecule, peptide_ksize,
             true_protein_coding_fasta):
    from khtools.partition_reads import cli

    runner = CliRunner()
    result = runner.invoke(cli,
                           ['--peptide-ksize', peptide_ksize,
                            "--peptides-are-bloom-filter",
                            '--molecule', molecule, peptide_bloom_filter,
                            reads])
    assert result.exit_code == 0
    assert result.output == true_protein_coding_fasta
