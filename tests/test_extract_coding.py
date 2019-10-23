from io import StringIO
import os
import warnings

from Bio.Seq import Seq
from click.testing import CliRunner
import pandas as pd
import pandas.util.testing as pdt
import pytest
import screed


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
    from khtools.extract_coding import three_frame_translation

    test = [str(x) for x in three_frame_translation(seq)]
    true = ['RLLNTDINNIRKIAI*L*ILFC', 'ACLILTSIILGKSQYNCKSCSV',
            'LA*Y*HQ*Y*ENRNITVNPVL']
    assert test == true


def test_three_frame_translation_no_stops(seq):
    from khtools.extract_coding import three_frame_translation_no_stops

    test = {k: str(v) for k, v in
            three_frame_translation_no_stops(seq).items()}
    true = {2: 'ACLILTSIILGKSQYNCKSCSV'}
    assert test == true


def test_six_frame_translation_no_stops(seq):
    from khtools.extract_coding import six_frame_translation_no_stops

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
def true_scores_path(data_folder, molecule, peptide_ksize):
    return os.path.join(data_folder, "extract_coding",
                            "SRR306838_GSM752691_hsa_br_F_1_trimmed_"
                            f"subsampled_n22__molecule-{molecule}_ksize-"
                            f"{peptide_ksize}.csv")


@pytest.fixture
def true_scores(true_scores_path):
    return pd.read_csv(true_scores_path, index_col=0)


@pytest.fixture
def true_protein_coding_fasta_path(data_folder):
    return os.path.join(data_folder, "extract_coding",
                        "true_protein_coding.fasta")

def test_score_reads(capsys, tmpdir, reads, peptide_bloom_filter, molecule,
                     peptide_ksize,  #true_scores,
                     true_scores_path,
                     true_protein_coding_fasta_path):
    from khtools.extract_coding import score_reads

    test = score_reads(reads, peptide_bloom_filter,
                       peptide_ksize=peptide_ksize, molecule=molecule)
    # pdt.assert_equal(test, true_scores)
    captured = capsys.readouterr()
    test_names = []
    for line in captured.out.splitlines():
        if line.startswith(">"):
            test_names.append(line.lstrip('>'))

    # Check that the proper sequences were output
    true_names = get_fasta_record_names(true_protein_coding_fasta_path)

    # Check that precision is high -- everything in "test" was truly coding
    assert all(test_name in true_names for test_name in test_names)

    # Check tqdm iterations
    assert '22it' in captured.err


def write_fasta_string_to_file(fasta_string, folder, prefix):
    test_fasta_filename = os.path.join(folder, prefix + '.fasta')
    with open(test_fasta_filename) as f:
        f.write(fasta_string)
    return test_fasta_filename


def get_fasta_record_names(fasta_path):
    names = []
    for record in screed.open(fasta_path):
        name = record['name']
        names.append(name)
    return set(names)


def test_cli_peptide_fasta(reads, peptide_fasta, molecule, peptide_ksize,
                           true_protein_coding_fasta_string):
    from khtools.extract_coding import cli

    runner = CliRunner()
    result = runner.invoke(cli,
                           ['--peptide-ksize', peptide_ksize,
                            '--molecule', molecule,
                            reads])
    assert result.exit_code == 0
    assert result.output == true_protein_coding_fasta_string


def test_cli_peptide_bloom_filter(reads, peptide_bloom_filter_path, molecule,
                                  peptide_ksize,
                                  true_protein_coding_fasta_string):
    from khtools.extract_coding import cli

    runner = CliRunner()
    result = runner.invoke(cli,
                           ['--peptide-ksize', peptide_ksize,
                            "--peptides-are-bloom-filter",
                            '--molecule', molecule, peptide_bloom_filter_path,
                            reads])
    assert result.exit_code == 0
    assert result.output == true_protein_coding_fasta_string
