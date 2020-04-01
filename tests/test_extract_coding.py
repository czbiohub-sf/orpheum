import os
import subprocess

import pandas as pd
import pandas.util.testing as pdt
import pytest
import screed

KHTOOLS = "khtools"
EXTRACT_CODING = "extract-coding"
CMD = KHTOOLS + " " + EXTRACT_CODING


@pytest.fixture
def low_complexity_seq():
    return "CCCCCCCCCACCACCACCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCA" \
           "CACACCCCCAACACCC"


@pytest.fixture(params=['seq', 'low_complexity_seq'])
def type_seq(request, seq, low_complexity_seq):
    if request.param == 'seq':
        return request.param, seq
    elif request.param == 'low_complexity_seq':
        return request.param, low_complexity_seq


def test_evaluate_is_fastp_low_complexity(type_seq):
    from khtools.extract_coding import evaluate_is_fastp_low_complexity

    seqtype, seq = type_seq

    test = evaluate_is_fastp_low_complexity(seq)
    if seqtype == 'seq':
        # regular sequence is not low complexity
        assert not test
    elif seqtype == 'low_complexity_seq':
        # low complexity sequence should evaluate to low complexity!
        assert test


@pytest.fixture
def empty_fasta(data_folder):
    return os.path.join(
        data_folder, 'empty_fasta.fasta')


@pytest.fixture
def true_scores_path(data_folder, alphabet, peptide_ksize):
    return os.path.join(
        data_folder, "extract_coding",
        "SRR306838_GSM752691_hsa_br_F_1_trimmed_"
        "subsampled_n22__alphabet-{}_ksize-".format(alphabet) +
        "{}.csv".format(peptide_ksize))


@pytest.fixture
def true_scores(true_scores_path):
    return pd.read_csv(true_scores_path)


@pytest.fixture
def true_protein_coding_fasta_path(data_folder):
    return os.path.join(data_folder, "extract_coding",
                        "true_protein_coding.fasta")


@pytest.fixture
def true_protein_coding_fasta_string(true_protein_coding_fasta_path):
    with open(true_protein_coding_fasta_path) as f:
        return f.read()


def get_fasta_record_names(fasta_path):
    names = []
    for record in screed.open(fasta_path):
        name = record['name']
        names.append(name)
    return set(names)


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
    print(test_scores)
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
