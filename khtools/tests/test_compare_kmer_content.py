"""
test_compare_kmer_content.py

Tests comparing k-mer content
"""
from io import StringIO

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

    test = kmer_comparison_table('seq1', nucleotide_seq1,
                                 'seq2', nucleotide_seq2,
                                 'nucleotide', ksizes=ksizes)
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
seq1,seq2,2,1.0,botvinnik
seq1,seq2,3,0.8888888888888888,botvinnik
seq1,seq2,4,0.75,botvinnik
seq1,seq2,2,1.0,dayhoff
seq1,seq2,3,0.8888888888888888,dayhoff
seq1,seq2,4,0.75,dayhoff
seq1,seq2,2,1.0,dayhoff_v2
seq1,seq2,3,0.8888888888888888,dayhoff_v2
seq1,seq2,4,0.75,dayhoff_v2
seq1,seq2,2,1.0,hydrophobic-polar
seq1,seq2,3,1.0,hydrophobic-polar
seq1,seq2,4,1.0,hydrophobic-polar
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
