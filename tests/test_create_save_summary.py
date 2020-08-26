import os
import warnings
import pytest
import pandas as pd
from sencha.create_save_summary import CreateSaveSummary
from sencha.constants_translate import (
    DEFAULT_JACCARD_THRESHOLD,
    LOW_COMPLEXITY_CATEGORIES,
)


@pytest.fixture()
def coding_scores_empty():
    return []


@pytest.fixture()
def coding_scores_nonempty():
    # Make fake dataframe with "read_id" and "category" columns only
    # for testing
    coding_scores = [
        ["read0", 0.9, 0, "Coding", 0, ""],
        ["read0", 0.9, 0, "Coding", 0, ""],
        ["read0", 0.9, 0, "Coding", 0, ""],
        ["read0", 0.9, 0, "Coding", 0, ""],
        ["read0", 0.9, 0, "Coding", 0, ""],
        ["read0", 0.9, 0, "Coding", 0, ""],
        ["read1", 0.9, 0, "Coding", 0, ""],
        ["read1", 0.9, 0, "Coding", 0, ""],
        ["read1", 0.9, 0, "Coding", 0, ""],
        ["read1", 0.9, 0, "Coding", 0, ""],
        ["read1", 0.9, 0, "Coding", 0, ""],
        ["read2", 0.9, 0, "Coding", 0, ""],
        ["read2", 0.9, 0, "Coding", 0, ""],
        ["read2", 0.9, 0, "Coding", 0, ""],
        ["read2", 0.9, 0, "Coding", 0, ""],
        ["read3", 0.9, 0, "Coding", 0, ""],
        ["read3", 0.9, 0, "Coding", 0, ""],
        ["read3", 0.9, 0, "Coding", 0, ""],
        ["read4", 0.9, 0, "Coding", 0, ""],
        ["read4", 0.9, 0, "Coding", 0, ""],
        ["read5", 0.9, 0, "Coding", 0, ""],
        ["read5", 0.9, 0, "Coding", 0, ""],
        ["read6", 0.9, 0, "Coding", 0, ""],
        ["read7", 0.9, 0, "Coding", 0, ""],
        ["read8", 0.9, 0, "Coding", 0, ""],
        ["read9", 0.9, 0, "Coding", 0, ""],
        ["read10", 0.9, 0, "Coding", 0, ""],
    ]
    return coding_scores


@pytest.fixture
def true_scores_path(data_folder):
    true_scores_path = os.path.join(
        data_folder,
        "translate",
        "SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_n22__"
        "alphabet-protein_ksize-7.csv",
    )
    return true_scores_path


@pytest.fixture
def true_scores_parquet(data_folder):
    true_scores_parquet = os.path.join(
        data_folder,
        "translate",
        "SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_n22__"
        "alphabet-protein_ksize-7.parquet",
    )
    return true_scores_parquet


@pytest.fixture
def single_alphabet_ksize_true_scores(true_scores_path):
    df = pd.read_csv(true_scores_path)
    return df.values.tolist()


def test_maybe_write_json_summary_empty(
    coding_scores_empty, alphabet, peptide_bloom_filter_path, peptide_ksize
):
    create_ss = CreateSaveSummary(
        ["nonexistent.fa"],
        True,
        True,
        True,
        peptide_bloom_filter_path,
        alphabet,
        peptide_ksize,
        DEFAULT_JACCARD_THRESHOLD,
        coding_scores_empty,
    )
    summary = create_ss.maybe_write_json_summary()
    assert summary["input_files"] == ["nonexistent.fa"]
    assert summary["jaccard_info"]["count"] == 0


def test_get_n_translated_frames_per_read(
    coding_scores_nonempty, alphabet, peptide_bloom_filter_path, peptide_ksize
):
    create_ss = CreateSaveSummary(
        ["nonexistent.fa"],
        True,
        True,
        True,
        peptide_bloom_filter_path,
        alphabet,
        peptide_ksize,
        DEFAULT_JACCARD_THRESHOLD,
        coding_scores_nonempty,
    )
    percentages, histogram = create_ss.get_n_translated_frames_per_read()
    assert histogram == {
        "Number of reads with 1 putative protein-coding translations": 5,
        "Number of reads with 2 putative protein-coding translations": 2,
        "Number of reads with 6 putative protein-coding translations": 1,
        "Number of reads with 5 putative protein-coding translations": 1,
        "Number of reads with 4 putative protein-coding translations": 1,
        "Number of reads with 3 putative protein-coding translations": 1,
    }
    assert percentages == {
        "Number of reads with 1 putative protein-coding translations": 45.45454545454545,
        "Number of reads with 2 putative protein-coding translations": 18.181818181818183,
        "Number of reads with 6 putative protein-coding translations": 9.090909090909092,
        "Number of reads with 5 putative protein-coding translations": 9.090909090909092,
        "Number of reads with 4 putative protein-coding translations": 9.090909090909092,
        "Number of reads with 3 putative protein-coding translations": 9.090909090909092,
    }


def test_get_n_per_coding_category(
    coding_scores_nonempty,
    alphabet,
    peptide_bloom_filter_path,
    peptide_ksize,
    jaccard_threshold,
):
    from sencha.sequence_encodings import ALIAS_TO_ALPHABET

    data = [
        ["read1", 0.9, 0, "Non-coding", 0, ""],
        ["read1", 0.9, 0, "Coding", 0, ""],
        ["read1", 0.9, 0, "Non-coding", 0, ""],
        ["read2", 0.9, 0, "Translation frame has stop codon(s)", 0, ""],
        ["read3", 0.9, 0, "Coding", 0, ""],
        ["read4", 0.9, 0, "Non-coding", 0, ""],
        ["read5", 0.9, 0, "Low complexity nucleotide", 0, ""],
        ["read6", 0.9, 0, "Read length was shorter than 3 * peptide k-mer size", 0, ""],
        ["read7", 0.9, 0, LOW_COMPLEXITY_CATEGORIES[alphabet], 0, ""],
    ]

    create_ss = CreateSaveSummary(
        ["nonexistent.fa"],
        True,
        True,
        True,
        peptide_bloom_filter_path,
        alphabet,
        peptide_ksize,
        jaccard_threshold,
        data,
    )

    test_counts, test_percentages = create_ss.get_n_per_coding_category()
    canonical_alphabet = ALIAS_TO_ALPHABET[alphabet]
    # read1 and read3 are coding, there is zero too_short_peptide
    true_counts = {
        "Translation is shorter than peptide k-mer size + 1": 0.0,
        "Translation frame has stop codon(s)": 14.285714285714286,
        "Coding": 28.571428571428573,
        "Non-coding": 14.285714285714286,
        "Low complexity nucleotide": 14.285714285714286,
        "Read length was shorter than 3 * peptide k-mer size": 14.285714285714286,
        f"Low complexity peptide in {canonical_alphabet} alphabet": 14.285714285714286,
    }
    true_percentages = {
        "Translation is shorter than peptide k-mer size + 1": 0,
        "Translation frame has stop codon(s)": 1,
        "Coding": 2,
        "Non-coding": 1,
        "Low complexity nucleotide": 1,
        "Read length was shorter than 3 * peptide k-mer size": 1,
        f"Low complexity peptide in {canonical_alphabet} alphabet": 1,
    }
    assert test_counts == true_counts
    assert test_percentages == true_percentages


def test_generate_coding_summary(reads, data_folder, single_alphabet_ksize_true_scores):
    create_ss = CreateSaveSummary(
        reads,
        True,
        True,
        True,
        "bloom_filter.nodegraph",
        "protein",
        7,
        0.5,
        single_alphabet_ksize_true_scores,
    )
    test_summary = create_ss.generate_coding_summary()
    true_summary = {
        "input_files": ["SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_n22.fq"],
        "jaccard_info": {
            "count": 44,
            "mean": 0.085830733808675,
            "std": 0.2503210514088884,
            "min": 0.0,
            "25%": 0.0,
            "50%": 0.0,
            "75%": 0.05882352941176471,
            "max": 1.0,
        },
        "categorization_counts": {
            "Translation is shorter than peptide k-mer size + 1": 0,
            "Translation frame has stop codon(s)": 3,
            "Coding": 3,
            "Non-coding": 14,
            "Low complexity nucleotide": 0,
            "Read length was shorter than 3 * peptide k-mer size": 2,
            "Low complexity peptide in protein20 alphabet": 1,
        },
        "categorization_percentages": {
            "Translation is shorter than peptide k-mer size + 1": 0.0,
            "Translation frame has stop codon(s)": 13.043478260869565,
            "Coding": 13.043478260869565,
            "Non-coding": 60.869565217391305,
            "Low complexity nucleotide": 0.0,
            "Read length was shorter than 3 * peptide k-mer size": 8.695652173913043,
            "Low complexity peptide in protein20 alphabet": 4.3478260869565215,
        },
        "histogram_n_coding_frames_per_read": {
            "Number of reads with 1 putative protein-coding translations": 3
        },
        "histogram_n_coding_frames_per_read_percentages": {
            "Number of reads with 1 putative protein-coding translations": 100.0
        },
        "peptide_bloom_filter": "bloom_filter.nodegraph",
        "peptide_alphabet": "protein",
        "peptide_ksize": 7,
        "jaccard_threshold": 0.5,
    }

    assert test_summary == true_summary


def test_maybe_write_csv(reads, single_alphabet_ksize_true_scores, true_scores_path):
    create_ss = CreateSaveSummary(
        reads,
        true_scores_path,
        True,
        True,
        "bloom_filter.nodegraph",
        "protein",
        7,
        0.5,
        single_alphabet_ksize_true_scores,
    )
    create_ss.maybe_write_csv()


def test_maybe_write_parquet(
    reads, single_alphabet_ksize_true_scores, true_scores_parquet
):
    create_ss = CreateSaveSummary(
        reads,
        True,
        true_scores_parquet,
        True,
        "bloom_filter.nodegraph",
        "protein",
        7,
        0.5,
        single_alphabet_ksize_true_scores,
    )
    create_ss.maybe_write_parquet()


def test_make_empty_coding_categories(single_alphabet_ksize_true_scores):
    create_ss = CreateSaveSummary(
        ["nonexistent.fa"],
        True,
        True,
        True,
        "bloom_filter.nodegraph",
        "protein",
        7,
        0.5,
        single_alphabet_ksize_true_scores,
    )
    test_coding_categories = {
        "Translation is shorter than peptide k-mer size + 1": 0,
        "Translation frame has stop codon(s)": 0,
        "Coding": 0,
        "Non-coding": 0,
        "Low complexity nucleotide": 0,
        "Read length was shorter than 3 * peptide k-mer size": 0,
        "Low complexity peptide in protein20 alphabet": 0,
    }
    true_coding_categories = create_ss.make_empty_coding_categories()
    assert true_coding_categories == test_coding_categories
