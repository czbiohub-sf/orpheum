import os
import sys

import sencha.fasta_utils as fasta_utils


def test_calculate_chunksize():
    tota_jobs_todo = 100
    processes = 1
    obtained = fasta_utils.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = fasta_utils.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained


def test_batch_iterator():
    description = "test"
    sequence = "seq"
    fasta_utils.batch_iterator(sys.stdout, description, sequence)


def test_split_fasta_files():
    description = "test"
    sequence = "seq"
    fasta_utils.split_fasta_files(sys.stdout, description, sequence)


def test_maybe_open_fastas():
    pass


def test_maybe_close_fastas():
    pass


def test_write_fasta(capsys):
    description = "test"
    sequence = "seq"
    fasta_utils.write_fasta(sys.stdout, description, sequence)
    captured = capsys.readouterr()
    assert captured.out == ">test\nseq\n"


def test_maybe_write_fasta(tmpdir, capsys):
    # Check if file handle is stdout
    description = "test"
    sequence = "seq"
    fasta_utils.maybe_write_fasta(sys.stdout, description, sequence)
    captured = capsys.readouterr()
    assert captured.out in ">test\nseq\n"
    # check if file handle is None
    fasta_utils.maybe_write_fasta(None, description, sequence)
    captured = capsys.readouterr()
    assert captured.out == ""
    fasta = os.path.join(tmpdir, "test_maybe_write_fasta.fasta")
    # check if file handle is a temporary fasta file
    fasta_utils.maybe_write_fasta(open(fasta, "w"), description, sequence)
    assert captured.out == ""


def test_open_and_announce(tmpdir, capsys):
    # Check if expected announcement is made
    fasta = os.path.join(tmpdir, "test_noncoding_nucleotide.fasta")
    fasta_utils.open_and_announce(fasta, "noncoding_nucleotide")
    captured = capsys.readouterr()
    expected = "Writing nucleotide sequence from reads WITHOUT matches to protein-coding peptides to {}\n".format(
        fasta
    )
    assert captured.out in expected
