import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import pytest

import khtools.extract_coding as ec
import khtools.constants_extract_coding as constants_ec
import khtools.constants_bloom_filter as constants_bf


KHTOOLS = "khtools"
EXTRACT_CODING = "extract-coding"
CMD = KHTOOLS + " " + EXTRACT_CODING


@pytest.fixture()
def extract_coding_class(tmpdir, reads, peptide_fasta):
    args = dict(
        reads=[reads],
        peptides=peptide_fasta,
        peptide_ksize=None,
        save_peptide_bloom_filter=True,
        peptides_are_bloom_filter=False,
        jaccard_threshold=None,
        alphabet='protein',
        csv=False,
        json_summary=False,
        coding_nucleotide_fasta=os.path.join(
            tmpdir, "coding_nucleotide_fasta.fa"),
        noncoding_nucleotide_fasta=os.path.join(
            tmpdir, "noncoding_nucleotide_fasta.fa"),
        low_complexity_nucleotide_fasta=os.path.join(
            tmpdir, "low_complexity_nucleotide_fasta.fa"),
        low_complexity_peptide_fasta=os.path.join(
            tmpdir, "low_complexity_peptide_fasta.fa"),
        tablesize=constants_bf.DEFAULT_MAX_TABLESIZE,
        n_tables=constants_bf.DEFAULT_N_TABLES,
        long_reads=False,
        verbose=True)
    ec_obj = ec.ExtractCoding(args)
    return ec_obj


@pytest.fixture
def translations_for_single_seq():
    return {
        2: 'ACLILTSIILGKSQYNCKSCSV',
        -2: 'TEQDLQLYCDFPNIIDVSIKQA',
        -3: 'QNRIYSYIAIFLILLMSVLSK'
    }


def test_get_jaccard_threshold():
    assert ec.get_jaccard_threshold(
        None, "") == constants_ec.DEFAULT_JACCARD_THRESHOLD
    assert ec.get_jaccard_threshold(0.5, "") == 0.5

    assert ec.get_jaccard_threshold(
        None, "protein") == constants_ec.DEFAULT_JACCARD_THRESHOLD
    assert ec.get_jaccard_threshold(
        None, "dayhoff") == constants_ec.DEFAULT_JACCARD_THRESHOLD
    assert ec.get_jaccard_threshold(
        None, "hp") == constants_ec.DEFAULT_HP_JACCARD_THRESHOLD


def test_evaluate_is_fastp_low_complexity(type_seq):
    seqtype, seq = type_seq
    test = ec.evaluate_is_fastp_low_complexity(seq)
    if seqtype == 'seq':
        # regular sequence is not low complexity
        assert not test
    elif seqtype == 'low_complexity_seq':
        # low complexity sequence should evaluate to low complexity!
        assert test


def test_compute_fastp_complexity(type_seq):
    seqtype, seq = type_seq
    test = ec.compute_fastp_complexity(seq)
    if seqtype == 'seq':
        # regular sequence is not low complexity
        np.testing.assert_almost_equal(test, 0.74, 0.001)
    elif seqtype == 'low_complexity_seq':
        # low complexity sequence should evaluate to low complexity!
        np.testing.assert_almost_equal(test, 0.26, 0.001)


def test_evaluate_is_kmer_low_complexity(type_seq):
    seqtype, seq = type_seq
    test = ec.evaluate_is_kmer_low_complexity(seq, 7)
    if seqtype == 'seq':
        # regular sequence is not low complexity
        assert test is False
    elif seqtype == 'low_complexity_seq':
        # low complexity sequence should evaluate to low complexity!
        assert test is True


def test_compute_kmer_complexity(type_seq):
    seqtype, seq = type_seq
    test = ec.compute_kmer_complexity(seq, 7)
    if seqtype == 'seq':
        # regular sequence is not low complexity
        assert test == 30.5
    elif seqtype == 'low_complexity_seq':
        # low complexity sequence should evaluate to low complexity!
        assert test == 35.0


def test_write_fasta(capsys):
    description = "test"
    sequence = "seq"
    ec.write_fasta(sys.stdout, description, sequence)
    captured = capsys.readouterr()
    assert captured.out == ">test\nseq\n"


def test_maybe_write_fasta(tmpdir, capsys, extract_coding_class):
    # Check if file handle is stdout
    description = "test"
    sequence = "seq"
    extract_coding_class.maybe_write_fasta(sys.stdout, description, sequence)
    captured = capsys.readouterr()
    assert captured.out == ">test\nseq\n"
    # check if file handle is None
    extract_coding_class.maybe_write_fasta(None, description, sequence)
    captured = capsys.readouterr()
    assert captured.out == ""
    fasta = os.path.join(tmpdir, 'test_maybe_write_fasta.fasta')
    # check if file handle is a temporary fasta file
    extract_coding_class.maybe_write_fasta(
        open(fasta, "w"), description, sequence)
    assert captured.out == ""


def test_open_and_announce(tmpdir, capsys, extract_coding_class):
    # Check if expected announcement is made
    fasta = os.path.join(tmpdir, 'test_noncoding_nucleotide.fasta')
    extract_coding_class.open_and_announce(fasta, "noncoding_nucleotide")
    captured = capsys.readouterr()
    expected = \
        'Writing nucleotide sequence from reads WITHOUT matches to protein-coding peptides to {}\n'.format(fasta)
    assert captured.out in expected


def test_maybe_open_fastas(tmpdir, capsys, extract_coding_class):
    # Check if file handle is stdout
    extract_coding_class.set_coding_scores_all_files()
    assert len(extract_coding_class.fastas) == 4
    assert len(extract_coding_class.file_handles) == 4
    fastas = [
        extract_coding_class.noncoding_nucleotide_fasta,
        extract_coding_class.coding_nucleotide_fasta,
        extract_coding_class.low_complexity_peptide_fasta,
        extract_coding_class.low_complexity_nucleotide_fasta]
    seqtypes = list(extract_coding_class.file_handles.keys())
    for key, value in extract_coding_class.fastas.items():
        assert value in fastas
        assert key in seqtypes


def test_ec_get_jaccard_threshold(extract_coding_class):
    assert extract_coding_class.get_jaccard_threshold(
    ) == constants_ec.DEFAULT_JACCARD_THRESHOLD


def test_score_single_translation(
        translations_for_single_seq, extract_coding_class):
    fraction_in_peptide_db, n_kmers = \
        extract_coding_class.score_single_translation(
            translations_for_single_seq)
    np.testing.assert_almost_equal(
        fraction_in_peptide_db, 0.19, 0.05)
    assert n_kmers == 82


def test_get_peptide_meta(
        extract_coding_class, translations_for_single_seq):

    # Convert to BioPython sequence object for translation
    result = \
        extract_coding_class.get_peptide_meta(
            translations_for_single_seq)
    expected = (
        {-3: 0.0, -2: 1.0, 2: 0.0},
        {-3: 15, -2: 16, 2: 16},
        {-3: False, -2: False, 2: False})
    assert result == expected


def test_get_coding_score_line(
        extract_coding_class, translations_for_single_seq):

    # Convert to BioPython sequence object for translation
    result = \
        extract_coding_class.get_coding_score_line(
            "description",
            0.5,
            40,
            "test special case")
    assert result == ['description', 0.5, 40, 'test special case']
    result = \
        extract_coding_class.get_coding_score_line(
            "description",
            1.0,
            40,
            "test special case")
    assert result == ['description', 1.0, 40, "test special case"]
    result = \
        extract_coding_class.get_coding_score_line(
            "description",
            1.0,
            40,
            None)
    assert result == ['description', 1.0, 40, 'Coding']


def test_cli_peptide_fasta(reads, peptide_fasta, alphabet, peptide_ksize,
                           true_protein_coding_fasta_string):
    proc = subprocess.Popen(
        CMD + " --peptide-ksize {}".format(peptide_ksize) +
        " --alphabet " + alphabet + " {} {}".format(peptide_fasta, reads),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    assert proc.returncode == 0
    # true string is contained in the output
    assert true_protein_coding_fasta_string in str(stdout, 'UTF-8')

    # Make sure "Writing extract_coding summary to" didn't get accidentally
    # written to stdout instead of stderr
    assert 'Writing extract_coding summary to' \
           not in true_protein_coding_fasta_string


def test_cli_bad_jaccard_threshold_float(reads, peptide_fasta):
    proc = subprocess.Popen(
        [KHTOOLS,
         EXTRACT_CODING,
         "--jaccard-threshold",
         "3.14",
         peptide_fasta, reads],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    assert proc.returncode == 2
    error_message = 'Error: Invalid value for "--jaccard-threshold": --jaccard-threshold needs to be a number between 0 and 1, but was provided'
    assert error_message in str(stdout, 'UTF-8')


def test_cli_bad_jaccard_threshold_string(reads, peptide_fasta):
    proc = subprocess.Popen(
        [KHTOOLS,
         EXTRACT_CODING,
         "--jaccard-threshold",
         "beyonce",
         peptide_fasta, reads],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 2
    error_message = 'Error: Invalid value for "--jaccard-threshold": beyonce is not a valid floating point value'
    assert error_message in str(stdout, 'UTF-8')


def test_cli_peptide_bloom_filter(reads, peptide_bloom_filter_path, alphabet,
                                  peptide_ksize,
                                  true_protein_coding_fasta_string):
    proc = subprocess.Popen(
        CMD + " --peptide-ksize" +
        " {} ".format(peptide_ksize) +
        "--peptides-are-bloom-filter --alphabet {} ".format(alphabet) +
        "{} {}".format(peptide_bloom_filter_path, reads),
        shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0
    assert true_protein_coding_fasta_string in str(stdout, 'UTF-8')


def test_cli_csv(tmpdir, reads, peptide_bloom_filter_path, alphabet,
                 peptide_ksize, true_protein_coding_fasta_string, true_scores):

    csv = os.path.join(tmpdir, 'coding_scores.csv')

    proc = subprocess.Popen(
        CMD +
        " --peptide-ksize {} ".format(peptide_ksize) +
        "--csv " + csv + " --peptides-are-bloom-filter " +
        "--alphabet " + alphabet + " {} {}".format(
            peptide_bloom_filter_path, reads),
        shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0
    assert true_protein_coding_fasta_string in str(stdout, 'UTF-8')
    assert os.path.exists(csv)

    # the CLI adds the filename to the scoring dataframe
    true = true_scores.copy()
    true['filename'] = reads

    test_scores = pd.read_csv(csv)
    pdt.assert_equal(test_scores, true)


def test_cli_coding_nucleotide_fasta(tmpdir, reads, peptide_fasta):
    coding_nucleotide_fasta = os.path.join(tmpdir, 'coding_nucleotides.fasta')

    result = subprocess.Popen(
        CMD +
        " --coding-nucleotide-fasta" +
        " {} {} {}".format(coding_nucleotide_fasta, peptide_fasta, reads),
        shell=True)
    stdout, stderr = result.communicate()
    assert result.returncode == 0
    assert os.path.exists(coding_nucleotide_fasta)


def test_cli_noncoding_fasta(tmpdir, reads, peptide_fasta):
    noncoding_nucleotide_fasta = os.path.join(tmpdir,
                                              'noncoding_nucleotides.fasta')

    result = subprocess.Popen(
        CMD +
        " --noncoding-nucleotide-fasta" +
        " {} {} {}".format(noncoding_nucleotide_fasta, peptide_fasta, reads),
        shell=True)
    stdout, stderr = result.communicate()
    assert result.returncode == 0
    assert os.path.exists(noncoding_nucleotide_fasta)


def test_cli_low_complexity_nucleotide(tmpdir, reads, peptide_fasta):
    low_complexity_nucleotide_fasta = os.path.join(
        tmpdir, 'low_complexity_nucleotide.fasta')

    result = subprocess.Popen(
        CMD +
        " --low-complexity-nucleotide-fasta" +
        " {} {} {}".format(
            low_complexity_nucleotide_fasta, peptide_fasta, reads),
        shell=True)
    stdout, stderr = result.communicate()
    assert result.returncode == 0
    assert os.path.exists(low_complexity_nucleotide_fasta)


def test_cli_low_complexity_peptide(
        tmpdir,
        reads,
        peptide_fasta):
    low_complexity_peptide_fasta = os.path.join(tmpdir,
                                                'low_complexity_peptide.fasta')

    result = subprocess.Popen(
        CMD + " --low-complexity-peptide-fasta" +
        " {} {} {}".format(low_complexity_peptide_fasta, peptide_fasta, reads),
        shell=True)
    stdout, stderr = result.communicate()
    assert result.returncode == 0
    assert os.path.exists(low_complexity_peptide_fasta)


def test_cli_json_summary(tmpdir, reads, peptide_fasta):
    json_summary = os.path.join(tmpdir, 'coding_summary.json')

    result = subprocess.Popen(
        CMD + " --json-summary" + " {} {} {}".format(
            json_summary, peptide_fasta, reads),
        shell=True)
    stdout, stderr = result.communicate()
    assert result.returncode == 0
    assert os.path.exists(json_summary)


def test_cli_empty_fasta_json_summary(tmpdir, empty_fasta, peptide_fasta):

    json_summary = os.path.join(tmpdir, 'coding_summary.json')

    result = subprocess.Popen(
        CMD + " --json-summary" + " {} {} {}".format(
            json_summary, peptide_fasta, empty_fasta),
        shell=True)
    stdout, stderr = result.communicate()
    assert result.returncode == 0
    assert os.path.exists(json_summary)
