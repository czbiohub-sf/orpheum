import os

import pandas as pd
import pytest

"""
conftest.py contains fixtures or functions-turned-variables that can be
used in any test
"""
from sencha.constants_index import (
    DEFAULT_PROTEIN_KSIZE,
    DEFAULT_DAYHOFF_KSIZE,
    DEFAULT_HP_KSIZE,
)


@pytest.fixture
def reads(data_folder):
    return os.path.join(
        data_folder, "SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_n22.fq"
    )


@pytest.fixture
def jaccard_threshold(alphabet):
    from sencha.translate import get_jaccard_threshold

    threshold = get_jaccard_threshold(None, alphabet)
    return threshold


@pytest.fixture
def seq():
    s = "CGCTTGCTTAATACTGACATCAATAATATTAGGAAAATCGCAATATAACTGTAAATCCTGTTCTGTC"
    return s


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), "./data")


@pytest.fixture
def peptide_fasta(data_folder):
    filename = os.path.join(
        data_folder, "index", "Homo_sapiens.GRCh38.pep.subset.fa.gz"
    )
    return filename


@pytest.fixture
def adversarial_peptide_fasta(data_folder):
    filename = os.path.join(
        data_folder, "index", "Homo_sapiens.GRCh38.pep.first1000lines.fa"
    )
    return filename


@pytest.fixture(params=["normal", "adversarial"])
def variable_peptide_fasta(request, peptide_fasta, adversarial_peptide_fasta):
    if request.param == "normal":
        return peptide_fasta
    else:
        return adversarial_peptide_fasta


@pytest.fixture
def peptides_dir(data_folder):
    foldername = os.path.join(data_folder, "index")
    return foldername


# Tie the alphabet name to its default ksize to make sure we keep getting the
# right sequences
@pytest.fixture(
    params=[
        ("protein", DEFAULT_PROTEIN_KSIZE),
        ("dayhoff", DEFAULT_DAYHOFF_KSIZE),
        # This k-mer size is too short for dayhoff and will
        # fail on the dataset
        pytest.param(("dayhoff", DEFAULT_PROTEIN_KSIZE), marks=pytest.mark.xfail),
        # Default peptide k-mer size for HP is 31 which
        # corresponds to a nucleotide read length of 90, and
        # is too long for this dataset
        pytest.param(("hydrophobic-polar", DEFAULT_HP_KSIZE), marks=pytest.mark.skip),
        # This is a wayyy too short k-mer size for the HP
        # alphabet and should definitely fail when scoring
        # reads
        pytest.param(
            ("hydrophobic-polar", DEFAULT_PROTEIN_KSIZE), marks=pytest.mark.xfail
        ),
    ],
    ids=[
        "protein_default_ksize",
        "dayhoff_default_ksize",
        "dayhoff_protein_ksize_xfail",
        "hp_default_ksize",
        "hp_protein_ksize_xfail",
    ],
)
def alphabet_ksize(request):
    return request.param


@pytest.fixture
def peptide_ksize(alphabet_ksize):
    return alphabet_ksize[1]


@pytest.fixture
def alphabet(alphabet_ksize):
    return alphabet_ksize[0]


@pytest.fixture
def peptide_bloom_filter_path(data_folder, alphabet, peptide_ksize):
    filename = os.path.join(
        data_folder,
        "index",
        "Homo_sapiens.GRCh38.pep.subset.alphabet-{}_".format(alphabet)
        + "ksize-{}.bloomfilter.nodegraph".format(peptide_ksize),
    )
    return filename


@pytest.fixture
def peptide_bloom_filter(
    peptide_bloom_filter_path, peptide_fasta, alphabet, peptide_ksize
):
    from sencha.index import load_nodegraph

    """Load bloom filter from path if exists, otherwise, make it"""
    try:
        return load_nodegraph(peptide_bloom_filter_path)
    except (OSError):
        from sencha.index import make_peptide_bloom_filter

        bloom_filter = make_peptide_bloom_filter(
            peptide_fasta, peptide_ksize, alphabet, tablesize=1e6
        )
        bloom_filter.save(peptide_bloom_filter_path)
        return bloom_filter


@pytest.fixture
def true_protein_coding_fasta_path(data_folder):
    return os.path.join(data_folder, "translate", "true_protein_coding.fasta")


@pytest.fixture
def true_protein_coding_fasta_string(true_protein_coding_fasta_path):
    with open(true_protein_coding_fasta_path) as f:
        return f.read()


@pytest.fixture
def low_complexity_seq():
    return (
        "CCCCCCCCCACCACCACCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCA"
        "CACACCCCCAACACCC"
    )


@pytest.fixture(params=["seq", "low_complexity_seq"])
def type_seq(request, seq, low_complexity_seq):
    if request.param == "seq":
        return request.param, seq
    elif request.param == "low_complexity_seq":
        return request.param, low_complexity_seq


@pytest.fixture
def empty_fasta(data_folder):
    return os.path.join(data_folder, "empty_fasta.fasta")


@pytest.fixture
def true_scores_path(data_folder, alphabet, peptide_ksize):
    return os.path.join(
        data_folder,
        "translate",
        "SRR306838_GSM752691_hsa_br_F_1_trimmed_"
        "subsampled_n22__alphabet-{}_ksize-".format(alphabet)
        + "{}.csv".format(peptide_ksize),
    )


@pytest.fixture
def true_scores(true_scores_path):
    return pd.read_csv(true_scores_path)
