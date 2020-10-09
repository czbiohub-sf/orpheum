import os
import sys

import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import pytest
from sencha.translate import cli
from click.testing import CliRunner

import sencha.translate as translate
import sencha.constants_translate as constants_translate
import sencha.constants_index as constants_index


KHTOOLS = "sencha"
TRANSLATE = "translate"
CMD = KHTOOLS + " " + TRANSLATE


@pytest.fixture()
def translate_class(tmpdir, reads, peptide_fasta):
    args = dict(
        reads=[reads],
        peptides=peptide_fasta,
        peptide_ksize=None,
        save_peptide_bloom_filter=True,
        peptides_are_bloom_filter=False,
        jaccard_threshold=None,
        alphabet="protein",
        csv=False,
        parquet=False,
        json_summary=False,
        coding_nucleotide_fasta=os.path.join(tmpdir, "coding_nucleotide_fasta.fa"),
        noncoding_nucleotide_fasta=os.path.join(
            tmpdir, "noncoding_nucleotide_fasta.fa"
        ),
        low_complexity_nucleotide_fasta=os.path.join(
            tmpdir, "low_complexity_nucleotide_fasta.fa"
        ),
        low_complexity_peptide_fasta=os.path.join(
            tmpdir, "low_complexity_peptide_fasta.fa"
        ),
        tablesize=constants_index.DEFAULT_MAX_TABLESIZE,
        n_tables=constants_index.DEFAULT_N_TABLES,
        long_reads=False,
        verbose=True,
    )
    translate_obj = translate.Translate(args)
    return translate_obj


@pytest.fixture
def translations_for_single_seq():
    return {
        2: "ACLILTSIILGKSQYNCKSCSV",
        -2: "TEQDLQLYCDFPNIIDVSIKQA",
        -3: "QNRIYSYIAIFLILLMSVLSK",
    }


def test_get_jaccard_threshold():
    assert (
        translate.get_jaccard_threshold(None, "")
        == constants_translate.DEFAULT_JACCARD_THRESHOLD
    )
    assert translate.get_jaccard_threshold(0.5, "") == 0.5

    assert (
        translate.get_jaccard_threshold(None, "protein")
        == constants_translate.DEFAULT_JACCARD_THRESHOLD
    )
    assert (
        translate.get_jaccard_threshold(None, "dayhoff")
        == constants_translate.DEFAULT_JACCARD_THRESHOLD
    )
    assert (
        translate.get_jaccard_threshold(None, "hp")
        == constants_translate.DEFAULT_HP_JACCARD_THRESHOLD
    )


def test_evaluate_is_fastp_low_complexity(type_seq):
    seqtype, seq = type_seq
    test = translate.evaluate_is_fastp_low_complexity(seq)
    if seqtype == "seq":
        # regular sequence is not low complexity
        assert not test
    elif seqtype == "low_complexity_seq":
        # low complexity sequence should evaluate to low complexity!
        assert test


def test_compute_fastp_complexity(type_seq):
    seqtype, seq = type_seq
    test = translate.compute_fastp_complexity(seq)
    if seqtype == "seq":
        # regular sequence is not low complexity
        np.testing.assert_almost_equal(test, 0.74, 0.001)
    elif seqtype == "low_complexity_seq":
        # low complexity sequence should evaluate to low complexity!
        np.testing.assert_almost_equal(test, 0.26, 0.001)


def test_evaluate_is_kmer_low_complexity(type_seq):
    seqtype, seq = type_seq
    test = translate.evaluate_is_kmer_low_complexity(seq, 7)
    if seqtype == "seq":
        # regular sequence is not low complexity
        assert test is False
    elif seqtype == "low_complexity_seq":
        # low complexity sequence should evaluate to low complexity!
        assert test is True


def test_compute_kmer_complexity(type_seq):
    seqtype, seq = type_seq
    test = translate.compute_kmer_complexity(seq, 7)
    if seqtype == "seq":
        # regular sequence is not low complexity
        assert test == 30.5
    elif seqtype == "low_complexity_seq":
        # low complexity sequence should evaluate to low complexity!
        assert test == 35.0


def test_write_fasta(capsys):
    description = "test"
    sequence = "seq"
    translate.write_fasta(sys.stdout, description, sequence)
    captured = capsys.readouterr()
    assert captured.out == ">test\nseq\n"


def test_maybe_write_fasta(tmpdir, capsys, translate_class):
    # Check if file handle is stdout
    description = "test"
    sequence = "seq"
    translate_class.maybe_write_fasta(sys.stdout, description, sequence)
    captured = capsys.readouterr()
    assert captured.out == ">test\nseq\n"
    # check if file handle is None
    translate_class.maybe_write_fasta(None, description, sequence)
    captured = capsys.readouterr()
    assert captured.out == ""
    fasta = os.path.join(tmpdir, "test_maybe_write_fasta.fasta")
    # check if file handle is a temporary fasta file
    translate_class.maybe_write_fasta(open(fasta, "w"), description, sequence)
    assert captured.out == ""


def test_open_and_announce(tmpdir, capsys, translate_class):
    # Check if expected announcement is made
    fasta = os.path.join(tmpdir, "test_noncoding_nucleotide.fasta")
    translate_class.open_and_announce(fasta, "noncoding_nucleotide")
    captured = capsys.readouterr()
    expected = "Writing nucleotide sequence from reads WITHOUT matches to protein-coding peptides to {}\n".format(
        fasta
    )
    assert captured.out in expected


def test_maybe_open_fastas(tmpdir, capsys, translate_class):
    # Check if file handle is stdout
    translate_class.set_coding_scores_all_files()
    assert len(translate_class.fastas) == 4
    assert len(translate_class.file_handles) == 4
    fastas = [
        translate_class.noncoding_nucleotide_fasta,
        translate_class.coding_nucleotide_fasta,
        translate_class.low_complexity_peptide_fasta,
        translate_class.low_complexity_nucleotide_fasta,
    ]
    seqtypes = list(translate_class.file_handles.keys())
    for key, value in translate_class.fastas.items():
        assert value in fastas
        assert key in seqtypes


def test_translate_get_jaccard_threshold(translate_class):
    assert (
        translate_class.get_jaccard_threshold()
        == constants_translate.DEFAULT_JACCARD_THRESHOLD
    )


def test_score_single_translation(translations_for_single_seq, translate_class):
    fraction_in_peptide_db, n_kmers = translate_class.score_single_translation(
        translations_for_single_seq
    )
    np.testing.assert_almost_equal(fraction_in_peptide_db, 0.19, 0.05)
    assert n_kmers == 82


def test_get_coding_score_line(translate_class, translations_for_single_seq):

    # Convert to BioPython sequence object for translation
    result = translate_class.get_coding_score_line(
        "description", 0.5, 40, "test special case", -2
    )
    assert result == ["description", 0.5, 40, "test special case", -2]
    result = translate_class.get_coding_score_line(
        "description", 1.0, 40, "test special case", 1
    )
    assert result == ["description", 1.0, 40, "test special case", 1]
    result = translate_class.get_coding_score_line("description", 1.0, 40, None, -3)
    assert result == ["description", 1.0, 40, "Coding", -3]


def test_cli_peptide_fasta(
    reads, peptide_fasta, alphabet, peptide_ksize, true_protein_coding_fasta_string
):

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--peptide-ksize",
            peptide_ksize,
            "--alphabet",
            alphabet,
            peptide_fasta,
            reads,
        ],
    )
    assert result.exit_code == 0
    # CliRunner jams together the stderr and stdout so just check if the
    # true string is contained in the output
    assert true_protein_coding_fasta_string in result.output

    # Make sure "Writing translate summary to" didn't get accidentally
    # written to stdout instead of stderr
    assert "Writing translate summary to" not in true_protein_coding_fasta_string


def test_cli_bad_jaccard_threshold_float(reads, peptide_fasta):

    runner = CliRunner()
    result = runner.invoke(cli, ["--jaccard-threshold", "3.14", peptide_fasta, reads])
    assert result.exit_code == 2
    error_message = "--jaccard-threshold needs to be a number between 0 and 1, but 3.14 was provided"
    assert error_message in result.output


def test_cli_bad_jaccard_threshold_string(reads, peptide_fasta):
    runner = CliRunner()
    result = runner.invoke(
        cli, ["--jaccard-threshold", "beyonce", peptide_fasta, reads]
    )
    assert result.exit_code == 2
    error_message = "beyonce is not a valid floating point value"
    assert error_message in result.output


def test_cli_peptide_bloom_filter(
    reads,
    peptide_bloom_filter_path,
    alphabet,
    peptide_ksize,
    true_protein_coding_fasta_string,
):
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--peptide-ksize",
            peptide_ksize,
            "--peptides-are-bloom-filter",
            "--alphabet",
            alphabet,
            peptide_bloom_filter_path,
            reads,
        ],
    )
    assert result.exit_code == 0
    assert true_protein_coding_fasta_string in result.output


def test_cli_csv(
    tmpdir,
    reads,
    peptide_bloom_filter_path,
    alphabet,
    peptide_ksize,
    true_protein_coding_fasta_string,
    true_scores,
):

    csv = os.path.join(tmpdir, "coding_scores.csv")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--peptide-ksize",
            peptide_ksize,
            "--csv",
            csv,
            "--peptides-are-bloom-filter",
            "--alphabet",
            alphabet,
            peptide_bloom_filter_path,
            reads,
        ],
    )
    assert result.exit_code == 0
    assert true_protein_coding_fasta_string in result.output
    assert os.path.exists(csv)

    # the CLI adds the filename to the scoring dataframe
    true = true_scores.copy()
    true["filename"] = reads

    test_scores = pd.read_csv(csv)
    pdt.assert_equal(test_scores, true)


def test_cli_parquet(
    tmpdir,
    reads,
    peptide_bloom_filter_path,
    alphabet,
    peptide_ksize,
    true_protein_coding_fasta_string,
    true_scores,
):

    parquet = os.path.join(tmpdir, "coding_scores.parquet")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--peptide-ksize",
            peptide_ksize,
            "--parquet",
            parquet,
            "--peptides-are-bloom-filter",
            "--alphabet",
            alphabet,
            peptide_bloom_filter_path,
            reads,
        ],
    )
    assert result.exit_code == 0
    assert true_protein_coding_fasta_string in result.output
    assert os.path.exists(parquet)

    # the CLI adds the filename to the scoring dataframe
    true = true_scores.copy()
    true["filename"] = reads

    test_scores = pd.read_parquet(parquet)
    pdt.assert_equal(test_scores, true)


def test_cli_coding_nucleotide_fasta(tmpdir, reads, peptide_fasta):

    coding_nucleotide_fasta = os.path.join(tmpdir, "coding_nucleotides.fasta")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["--coding-nucleotide-fasta", coding_nucleotide_fasta, peptide_fasta, reads],
    )
    assert result.exit_code == 0
    assert os.path.exists(coding_nucleotide_fasta)


def test_cli_noncoding_fasta(tmpdir, reads, peptide_fasta):

    noncoding_nucleotide_fasta = os.path.join(tmpdir, "noncoding_nucleotides.fasta")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--noncoding-nucleotide-fasta",
            noncoding_nucleotide_fasta,
            peptide_fasta,
            reads,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(noncoding_nucleotide_fasta)


def test_cli_low_complexity_nucleotide(tmpdir, reads, peptide_fasta):

    low_complexity_nucleotide_fasta = os.path.join(
        tmpdir, "low_complexity_nucleotide.fasta"
    )

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--low-complexity-nucleotide-fasta",
            low_complexity_nucleotide_fasta,
            peptide_fasta,
            reads,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(low_complexity_nucleotide_fasta)


def test_cli_low_complexity_peptide(tmpdir, reads, peptide_fasta):

    low_complexity_peptide_fasta = os.path.join(tmpdir, "low_complexity_peptide.fasta")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--low-complexity-peptide-fasta",
            low_complexity_peptide_fasta,
            peptide_fasta,
            reads,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(low_complexity_peptide_fasta)


def test_cli_json_summary(tmpdir, reads, peptide_fasta):

    json_summary = os.path.join(tmpdir, "coding_summary.json")

    runner = CliRunner()
    result = runner.invoke(cli, ["--json-summary", json_summary, peptide_fasta, reads])
    assert result.exit_code == 0
    assert os.path.exists(json_summary)


def test_cli_empty_fasta_json_summary(tmpdir, empty_fasta, peptide_fasta):

    json_summary = os.path.join(tmpdir, "coding_summary.json")

    runner = CliRunner()
    result = runner.invoke(
        cli, ["--json-summary", json_summary, peptide_fasta, empty_fasta]
    )
    assert result.exit_code == 0
    assert os.path.exists(json_summary)
