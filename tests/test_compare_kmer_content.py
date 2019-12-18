"""
test_compare_kmer_content.py

Tests comparing k-mer content
"""
from io import StringIO
import os

from click.testing import CliRunner
import pandas as pd
import pandas.util.testing as pdt
import pytest


# Fixtures are functions-turned-variables that can be used across multiple
# tests. conftest.py contains fixtures that can be used by any test file
@pytest.fixture
def folder():
    return "test-folder"


@pytest.fixture
def ksize():
    return 3


@pytest.fixture
def ksizes():
    return [2, 3, 4]


@pytest.fixture
def nucleotide_seq1():
    return "GATTACA"


@pytest.fixture
def nucleotide_seq2():
    return "GATTTTAAAACA"


@pytest.fixture
def peptide_seq1():
    return "SASHAFIERCE"


@pytest.fixture
def peptide_seq2():
    return "YASSSSASHAAAFIERCE"


def test_kmerize(nucleotide_seq1, ksize):
    from khtools.compare_kmer_content import kmerize

    test = kmerize(nucleotide_seq1, ksize)
    true = {'ACA', 'ATT', 'GAT', 'TAC', 'TTA'}
    assert test == true


def test_jaccardize(nucleotide_seq1):
    from khtools.compare_kmer_content import jaccardize

    test = jaccardize({1, 2, 3, 4}, {3, 4, 5, 6})
    true = 0.5
    assert test == true


def test_kmer_comparison_table(nucleotide_seq1, nucleotide_seq2, ksizes):
    from khtools.compare_kmer_content import kmer_comparison_table

    test = kmer_comparison_table('seq1',
                                 nucleotide_seq1,
                                 'seq2',
                                 nucleotide_seq2,
                                 'nucleotide',
                                 ksizes=ksizes)
    s = """id1,id2,ksize,jaccard,molecule
seq1,seq2,2,1.0,nucleotide
seq1,seq2,3,0.8,nucleotide
seq1,seq2,4,0.25,nucleotide
"""
    true = pd.read_csv(StringIO(s))
    pdt.assert_equal(test, true)


def test_compare_peptide_seqs(peptide_seq1, peptide_seq2, ksizes):
    from khtools.compare_kmer_content import compare_peptide_seqs
    id_seq1 = "seq1", peptide_seq1
    id_seq2 = "seq2", peptide_seq2

    test = compare_peptide_seqs(id_seq1, id_seq2, ksizes)

    s = """id1,id2,ksize,jaccard,molecule
seq1,seq2,2,1.0,protein
seq1,seq2,3,0.8888888888888888,protein
seq1,seq2,4,0.75,protein
seq1,seq2,2,1.0,dayhoff
seq1,seq2,3,0.8888888888888888,dayhoff
seq1,seq2,4,0.75,dayhoff
seq1,seq2,2,1.0,hp
seq1,seq2,3,1.0,hp
seq1,seq2,4,1.0,hp
seq1,seq2,2,1.0,botvinnik
seq1,seq2,3,0.8888888888888888,botvinnik
seq1,seq2,4,0.75,botvinnik
seq1,seq2,2,1.0,aa9
seq1,seq2,3,0.8888888888888888,aa9
seq1,seq2,4,0.75,aa9
seq1,seq2,2,1.0,gbmr4
seq1,seq2,3,1.0,gbmr4
seq1,seq2,4,0.7142857142857143,gbmr4
seq1,seq2,2,1.0,sdm12
seq1,seq2,3,0.8888888888888888,sdm12
seq1,seq2,4,0.75,sdm12
seq1,seq2,2,1.0,hsdm17
seq1,seq2,3,0.8888888888888888,hsdm17
seq1,seq2,4,0.75,hsdm17
"""
    true = pd.read_csv(StringIO(s))
    pdt.assert_equal(test, true)


def test_compare_nucleotide_seqs(nucleotide_seq1, nucleotide_seq2, ksizes):
    from khtools.compare_kmer_content import compare_nucleotide_seqs
    id_seq1 = "seq1", nucleotide_seq1
    id_seq2 = "seq2", nucleotide_seq2

    test = compare_nucleotide_seqs(id_seq1, id_seq2, ksizes)

    s = """id1,id2,ksize,jaccard,molecule
seq1,seq2,2,1.0,purine_pyrimidine
seq1,seq2,3,0.8,purine_pyrimidine
seq1,seq2,4,0.25,purine_pyrimidine
seq1,seq2,2,1.0,nucleotide
seq1,seq2,3,0.8,nucleotide
seq1,seq2,4,0.25,nucleotide
seq1,seq2,2,1.0,weak_strong
seq1,seq2,3,1.0,weak_strong
seq1,seq2,4,1.0,weak_strong
seq1,seq2,2,1.0,amino_keto
seq1,seq2,3,1.0,amino_keto
seq1,seq2,4,0.75,amino_keto
"""
    true = pd.read_csv(StringIO(s))
    pdt.assert_equal(test, true)


@pytest.fixture
def ksize_args():
    return ['--ksize-max', '3']


@pytest.fixture
def true_comparison_csv_kmax3():
    return '''id1,id2,ksize,jaccard,molecule
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,2,0.0,protein
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,3,0.0,protein
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,2,0.2631578947368421,botvinnik
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,3,0.0,botvinnik
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,2,0.4666666666666667,dayhoff
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,3,0.05263157894736842,dayhoff
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,2,0.4375,dayhoff_v2
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,3,0.0,dayhoff_v2
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,2,1.0,hydrophobic-polar
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,3,1.0,hydrophobic-polar
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,2,0.047619047619047616,protein
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,3,0.0,protein
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,2,0.2631578947368421,botvinnik
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,3,0.05,botvinnik
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,2,0.6428571428571429,dayhoff
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,3,0.2631578947368421,dayhoff
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,2,0.5,dayhoff_v2
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,3,0.21052631578947367,dayhoff_v2
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,2,1.0,hydrophobic-polar
SRR306838.10559374 Ibis_Run100924_C3PO:6:51:17601:17119/1 translation_frame: -2 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,3,1.0,hydrophobic-polar
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,2,0.14285714285714285,protein
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,3,0.0,protein
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,2,0.3157894736842105,botvinnik
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,3,0.0,botvinnik
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,2,0.7857142857142857,dayhoff
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,3,0.2631578947368421,dayhoff
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,2,0.625,dayhoff_v2
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,3,0.05263157894736842,dayhoff_v2
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,2,1.0,hydrophobic-polar
SRR306838.2740879 Ibis_Run100924_C3PO:6:13:11155:5248/1 translation_frame: -1 jaccard: 1.0,SRR306838.4880582 Ibis_Run100924_C3PO:6:23:17413:5436/1 translation_frame: 2 jaccard: 1.0,3,1.0,hydrophobic-polar'''


@pytest.fixture
def true_comparison_df(true_comparison_csv_kmax3):
    return pd.read_csv(StringIO(true_comparison_csv_kmax3))


def test_cli(true_protein_coding_fasta_path, ksize_args,
             true_comparison_csv_kmax3):
    from khtools.compare_kmer_content import cli

    runner = CliRunner()
    args = ksize_args + [true_protein_coding_fasta_path]
    result = runner.invoke(cli, args)
    assert result.exit_code == 0
    assert true_comparison_csv_kmax3 in result.output


def test_cli_no_input():
    from khtools.compare_kmer_content import cli

    runner = CliRunner()
    result = runner.invoke(cli)
    assert result.exit_code != 0
    error_message = "No sequence files provided! Argument 'fastas' is " \
                    "required!"
    assert error_message in result.exception.args[0]


def test_cli_parquet_no_csv(tmpdir, true_protein_coding_fasta_path, ksize_args,
                            true_comparison_csv_kmax3, true_comparison_df):
    from khtools.compare_kmer_content import cli

    parquet = os.path.join(tmpdir, 'coding_scores.parquet')

    runner = CliRunner()
    args = ksize_args + ['--parquet', parquet, '--no-csv',
                         true_protein_coding_fasta_path]
    result = runner.invoke(cli, args)
    assert result.exit_code == 0

    # Check that csv was not output to stdout
    assert true_comparison_csv_kmax3 not in result.output

    # Check existence of parquet file
    assert os.path.exists(parquet)

    # the CLI adds the filename to the scoring dataframe
    true = true_comparison_df.copy()

    test = pd.read_parquet(parquet)
    pdt.assert_equal(test, true)
