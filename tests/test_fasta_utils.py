import os
import sys

import sencha.fasta_utils as fasta_utils


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


def test_maybe_open_fastas(tmpdir, capsys):
    # Check if file handle is stdout
    fasta_utils.set_coding_scores_all_files()
    assert len(fasta_utils.fastas) == 4
    assert len(fasta_utils.file_handles) == 4
    fastas = [
        fasta_utils.noncoding_nucleotide_fasta,
        fasta_utils.coding_nucleotide_fasta,
        fasta_utils.low_complexity_peptide_fasta,
        fasta_utils.low_complexity_nucleotide_fasta,
    ]
    seqtypes = list(fasta_utils.file_handles.keys())
    for key, value in fasta_utils.fastas.items():
        assert value in fastas
        assert key in seqtypes
