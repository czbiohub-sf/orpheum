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
    from sencha.compare_kmer_content import kmerize

    test = kmerize(nucleotide_seq1, ksize)
    true = {"ACA", "ATT", "GAT", "TAC", "TTA"}
    assert test == true


def test_jaccardize(nucleotide_seq1):
    from sencha.compare_kmer_content import jaccardize

    test = jaccardize({1, 2, 3, 4}, {3, 4, 5, 6})
    true = 0.5
    assert test == true


def test_kmer_comparison_table(nucleotide_seq1, nucleotide_seq2, ksizes):
    from sencha.compare_kmer_content import kmer_comparison_table

    test = kmer_comparison_table(
        "seq1", nucleotide_seq1, "seq2", nucleotide_seq2, "nucleotide", ksizes=ksizes
    )
    s = """id1,id2,ksize,jaccard,alphabet
seq1,seq2,2,1.0,nucleotide
seq1,seq2,3,0.8,nucleotide
seq1,seq2,4,0.25,nucleotide
"""
    true = pd.read_csv(StringIO(s))
    pdt.assert_equal(test, true)


def test_compare_peptide_seqs(peptide_seq1, peptide_seq2, ksizes):
    from sencha.compare_kmer_content import compare_peptide_seqs

    id_seq1 = "seq1", peptide_seq1
    id_seq2 = "seq2", peptide_seq2

    test = compare_peptide_seqs(id_seq1, id_seq2, ksizes)
    data = [
        ["seq1", "seq2", 2, 1.0, "peptide20"],
        ["seq1", "seq2", 3, 0.8888888888888888, "peptide20"],
        ["seq1", "seq2", 4, 0.75, "peptide20"],
        ["seq1", "seq2", 2, 1.0, "hsdm17"],
        ["seq1", "seq2", 3, 0.8888888888888888, "hsdm17"],
        ["seq1", "seq2", 4, 0.75, "hsdm17"],
        ["seq1", "seq2", 2, 1.0, "sdm12"],
        ["seq1", "seq2", 3, 0.8888888888888888, "sdm12"],
        ["seq1", "seq2", 4, 0.75, "sdm12"],
        ["seq1", "seq2", 2, 1.0, "aa9"],
        ["seq1", "seq2", 3, 0.8888888888888888, "aa9"],
        ["seq1", "seq2", 4, 0.75, "aa9"],
        ["seq1", "seq2", 2, 1.0, "botvinnik8"],
        ["seq1", "seq2", 3, 0.8888888888888888, "botvinnik8"],
        ["seq1", "seq2", 4, 0.75, "botvinnik8"],
        ["seq1", "seq2", 2, 1.0, "dayhoff6"],
        ["seq1", "seq2", 3, 0.8888888888888888, "dayhoff6"],
        ["seq1", "seq2", 4, 0.75, "dayhoff6"],
        ["seq1", "seq2", 2, 1.0, "gbmr4"],
        ["seq1", "seq2", 3, 1.0, "gbmr4"],
        ["seq1", "seq2", 4, 0.7142857142857143, "gbmr4"],
        ["seq1", "seq2", 2, 1.0, "hp2"],
        ["seq1", "seq2", 3, 1.0, "hp2"],
        ["seq1", "seq2", 4, 1.0, "hp2"],
    ]

    true = pd.DataFrame(data, columns=["id1", "id2", "ksize", "jaccard", "alphabet"])
    pdt.assert_equal(test, true)


def test_compare_nucleotide_seqs(nucleotide_seq1, nucleotide_seq2, ksizes):
    from sencha.compare_kmer_content import compare_nucleotide_seqs

    id_seq1 = "seq1", nucleotide_seq1
    id_seq2 = "seq2", nucleotide_seq2

    test = compare_nucleotide_seqs(id_seq1, id_seq2, ksizes)

    s = """id1,id2,ksize,jaccard,alphabet
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
    return ["--ksize-max", "3"]


@pytest.fixture
def true_comparison_csv_kmax3():
    return """id1,id2,ksize,jaccard,alphabet
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,1.0,peptide20
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,1.0,peptide20
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,1.0,hsdm17
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,1.0,hsdm17
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,1.0,sdm12
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,1.0,sdm12
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,1.0,aa9
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,1.0,aa9
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,1.0,botvinnik8
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,1.0,botvinnik8
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,1.0,dayhoff6
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,1.0,dayhoff6
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,1.0,gbmr4
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,1.0,gbmr4
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,1.0,hp2
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,1.0,hp2
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.0,peptide20
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.0,peptide20
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.05263157894736842,hsdm17
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.0,hsdm17
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.17647058823529413,sdm12
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.0,sdm12
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.4666666666666667,aa9
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.0,aa9
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.2631578947368421,botvinnik8
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.0,botvinnik8
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.4666666666666667,dayhoff6
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.05263157894736842,dayhoff6
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.8,gbmr4
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.75,gbmr4
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,1.0,hp2
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,1.0,hp2
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.047619047619047616,peptide20
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.0,peptide20
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.15789473684210525,hsdm17
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.0,hsdm17
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.375,sdm12
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.0,sdm12
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.5333333333333333,aa9
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.05555555555555555,aa9
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.2631578947368421,botvinnik8
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.05,botvinnik8
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.6428571428571429,dayhoff6
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.2631578947368421,dayhoff6
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.6666666666666666,gbmr4
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.6,gbmr4
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,1.0,hp2
SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,1.0,hp2
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.0,peptide20
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.0,peptide20
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.05263157894736842,hsdm17
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.0,hsdm17
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.17647058823529413,sdm12
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.0,sdm12
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.4666666666666667,aa9
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.0,aa9
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.2631578947368421,botvinnik8
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.0,botvinnik8
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.4666666666666667,dayhoff6
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.05263157894736842,dayhoff6
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.8,gbmr4
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.75,gbmr4
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,1.0,hp2
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,1.0,hp2
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,1.0,peptide20
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,1.0,peptide20
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,1.0,hsdm17
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,1.0,hsdm17
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,1.0,sdm12
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,1.0,sdm12
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,1.0,aa9
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,1.0,aa9
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,1.0,botvinnik8
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,1.0,botvinnik8
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,1.0,dayhoff6
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,1.0,dayhoff6
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,1.0,gbmr4
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,1.0,gbmr4
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,1.0,hp2
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,1.0,hp2
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.14285714285714285,peptide20
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.0,peptide20
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.23809523809523808,hsdm17
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.05,hsdm17
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.5625,sdm12
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.2777777777777778,sdm12
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.5294117647058824,aa9
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.1,aa9
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.3157894736842105,botvinnik8
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.0,botvinnik8
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.7857142857142857,dayhoff6
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.2631578947368421,dayhoff6
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,0.8,gbmr4
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,0.5,gbmr4
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,1.0,hp2
SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,1.0,hp2
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.047619047619047616,peptide20
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.0,peptide20
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.15789473684210525,hsdm17
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.0,hsdm17
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.375,sdm12
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.0,sdm12
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.5333333333333333,aa9
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.05555555555555555,aa9
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.2631578947368421,botvinnik8
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.05,botvinnik8
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.6428571428571429,dayhoff6
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.2631578947368421,dayhoff6
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,0.6666666666666666,gbmr4
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,0.6,gbmr4
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,2,1.0,hp2
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.10559374__Ibis_Run100924_C3PO:6:51:17601:17119/1__translation_frame:-2__jaccard:1.0,3,1.0,hp2
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.14285714285714285,peptide20
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.0,peptide20
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.23809523809523808,hsdm17
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.05,hsdm17
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.5625,sdm12
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.2777777777777778,sdm12
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.5294117647058824,aa9
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.1,aa9
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.3157894736842105,botvinnik8
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.0,botvinnik8
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.7857142857142857,dayhoff6
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.2631578947368421,dayhoff6
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,0.8,gbmr4
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,0.5,gbmr4
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,2,1.0,hp2
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.2740879__Ibis_Run100924_C3PO:6:13:11155:5248/1__translation_frame:-1__jaccard:1.0,3,1.0,hp2
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,1.0,peptide20
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,1.0,peptide20
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,1.0,hsdm17
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,1.0,hsdm17
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,1.0,sdm12
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,1.0,sdm12
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,1.0,aa9
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,1.0,aa9
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,1.0,botvinnik8
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,1.0,botvinnik8
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,1.0,dayhoff6
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,1.0,dayhoff6
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,1.0,gbmr4
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,1.0,gbmr4
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,2,1.0,hp2
SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,SRR306838.4880582__Ibis_Run100924_C3PO:6:23:17413:5436/1__translation_frame:2__jaccard:1.0,3,1.0,hp2
"""


@pytest.fixture
def true_comparison_df(true_comparison_csv_kmax3):
    return pd.read_csv(StringIO(true_comparison_csv_kmax3))


def test_cli(true_protein_coding_fasta_path, ksize_args, true_comparison_csv_kmax3):
    from sencha.compare_kmer_content import cli

    runner = CliRunner()
    args = ksize_args + [true_protein_coding_fasta_path]
    result = runner.invoke(cli, args)
    assert result.exit_code == 0
    assert true_comparison_csv_kmax3 in result.output


def test_cli_no_input():
    from sencha.compare_kmer_content import cli

    runner = CliRunner()
    result = runner.invoke(cli)
    assert result.exit_code != 0
    error_message = "No sequence files provided! Argument 'fastas' is " "required!"
    assert error_message in result.exception.args[0]


def test_cli_parquet_no_csv(
    tmpdir,
    true_protein_coding_fasta_path,
    ksize_args,
    true_comparison_csv_kmax3,
    true_comparison_df,
):
    from sencha.compare_kmer_content import cli

    parquet = os.path.join(tmpdir, "coding_scores.parquet")

    runner = CliRunner()
    args = ksize_args + [
        "--parquet",
        parquet,
        "--no-csv",
        true_protein_coding_fasta_path,
    ]
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
