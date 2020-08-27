import os
import sys
import tempfile

import screed

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
    expected_batches = {0: [0, 1, 2], 1: [3, 4, 5], 2: [6, 7, 8], 3: [9, 10]}
    for i, batch in enumerate(fasta_utils.batch_iterator(iter(range(11)), 3)):
        assert batch == expected_batches[i]


def test_split_fasta_files(reads):
    temp_folder = tempfile.mkdtemp()
    filenames = fasta_utils.split_fasta_files(reads, 4, temp_folder)
    assert len(filenames) == 6
    expected_record_counts = {0: 4, 1: 4, 2: 4, 3: 4, 4: 4, 5: 3}
    for i, filename in enumerate(filenames):
        count = 0
        for record in screed.open(filename):
            count += 1
        assert count == expected_record_counts[i]


def test_maybe_open_fastas(tmpdir, capsys):
    description = "test"
    sequence = "seq"
    fasta = tmpdir + "test.fasta"
    fasta_utils.write_fasta(fasta, description, sequence)
    description = "test2"
    sequence = "seq2"
    fasta2 = tmpdir + "test2.fasta"
    fasta_utils.write_fasta(fasta, description, sequence)
    test_fastas = {"test_fasta": fasta, "test1_fasta": fasta2}
    announcement_dict = {"test_fasta": "Writing test1", "test1_fasta": "Writing test2"}
    file_handles = fasta_utils.maybe_open_fastas(test_fastas, announcement_dict)
    captured = capsys.readouterr()
    fasta_keys = list(test_fastas.keys())
    for key, value in file_handles.items():
        assert key in fasta_keys
        assert captured.out in announcement_dict[key]


def test_maybe_close_fastas(tmpdir):
    description = "test"
    sequence = "seq"
    fasta = tmpdir + "test.fasta"
    fasta_utils.write_fasta(fasta, description, sequence)
    description = "test2"
    sequence = "seq2"
    fasta2 = tmpdir + "test2.fasta"
    fasta_utils.write_fasta(fasta, description, sequence)
    test_fastas = {"test_fasta": fasta, "test1_fasta": fasta2}
    announcement_dict = {"test_fasta": "Writing test1", "test1_fasta": "Writing test2"}
    file_handles = fasta_utils.maybe_open_fastas(test_fastas, announcement_dict)
    fasta_utils.maybe_close_fastas(file_handles)


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
