import os
import warnings
import pytest
import pandas as pd
from khtools.assemble_coding_summary import AssembleSaveSummary
from khtools.constants_extract_coding import (
    DEFAULT_JACCARD_THRESHOLD, LOW_COMPLEXITY_CATEGORIES)


@pytest.fixture()
def coding_scores_empty():
    coding_scores = pd.DataFrame(
        columns=['read_id', 'jaccard_in_peptide_db',
                 'n_kmers', 'classification',
                 'filename'])
    return coding_scores


@pytest.fixture()
def coding_scores_nonempty():
    # Make fake dataframe with "read_id" and "classification" columns only
    # for testing
    data = {
        'read_id':
            ['read0', 'read0', 'read0', 'read0', 'read0', 'read0',
             'read1', 'read1', 'read1', 'read1', 'read1',
             'read2', 'read2', 'read2', 'read2',
             'read3', 'read3', 'read3',
             'read4', 'read4',
             'read5', 'read5',
             'read6', 'read7', 'read8', 'read9', 'read10'],
    }
    df = pd.DataFrame(data)
    df['classification'] = "Coding"
    return df


@pytest.fixture
def true_scores_path(data_folder):
    true_scores_path = os.path.join(
        data_folder, "extract_coding",
        "SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_n22__"
        "alphabet-protein_ksize-7.csv")
    return true_scores_path


@pytest.fixture
def single_alphabet_ksize_true_scores(true_scores_path):
    return pd.read_csv(true_scores_path)


def test_maybe_write_json_summary_empty(
        coding_scores_empty, alphabet,
        peptide_bloom_filter_path, peptide_ksize):
    assemble_ss = AssembleSaveSummary(
        ['nonexistent.fa'], True, True,
        peptide_bloom_filter_path,
        alphabet, peptide_ksize, DEFAULT_JACCARD_THRESHOLD)
    summary = assemble_ss.maybe_write_json_summary(coding_scores_empty)
    assert summary['input_files'] == ['nonexistent.fa']
    assert summary['jaccard_info']['count'] == 0


def test_get_n_translated_frames_per_read(
        coding_scores_nonempty, alphabet,
        peptide_bloom_filter_path, peptide_ksize):
    assemble_ss = AssembleSaveSummary(
        ['nonexistent.fa'], True, True,
        peptide_bloom_filter_path,
        alphabet, peptide_ksize, DEFAULT_JACCARD_THRESHOLD)
    percentages, histogram = assemble_ss.get_n_translated_frames_per_read(
        coding_scores_nonempty)
    assert histogram == {
        'Number of reads with 1 putative protein-coding translations': 5,
        'Number of reads with 2 putative protein-coding translations': 2,
        'Number of reads with 6 putative protein-coding translations': 1,
        'Number of reads with 5 putative protein-coding translations': 1,
        'Number of reads with 4 putative protein-coding translations': 1,
        'Number of reads with 3 putative protein-coding translations': 1}
    assert percentages == {
        'Number of reads with 1 putative protein-coding translations': 45.45454545454545,
        'Number of reads with 2 putative protein-coding translations': 18.181818181818183,
        'Number of reads with 6 putative protein-coding translations': 9.090909090909092,
        'Number of reads with 5 putative protein-coding translations': 9.090909090909092,
        'Number of reads with 4 putative protein-coding translations': 9.090909090909092,
        'Number of reads with 3 putative protein-coding translations': 9.090909090909092}


def test_get_n_per_coding_classification(
        coding_scores_nonempty, alphabet,
        peptide_bloom_filter_path, peptide_ksize, jaccard_threshold):
    from khtools.sequence_encodings import ALIAS_TO_ALPHABET
    assemble_ss = AssembleSaveSummary(
        ['nonexistent.fa'], True, True,
        peptide_bloom_filter_path,
        alphabet, peptide_ksize, jaccard_threshold)
    data = [
        ['read1', 'All translations shorter than peptide k-mer size + 1'],
        ['read2', 'All translation frames have stop codons'],
        ['read3', 'Coding'],
        ['read4', 'Non-coding'],
        ['read5', 'Low complexity nucleotide'],
        ['read6', 'Read length was shorter than 3 * peptide k-mer size'],
        ['read7', LOW_COMPLEXITY_CATEGORIES[alphabet]],
    ]
    df = pd.DataFrame(data, columns=['read_id', 'classification'])

    test_counts, test_percentages = \
        assemble_ss.get_n_per_coding_classification(df)
    canonical_alphabet = ALIAS_TO_ALPHABET[alphabet]
    true_counts = {
        'All translations shorter than peptide k-mer size + 1': 14.285714285714286,
        'All translation frames have stop codons': 14.285714285714286,
        'Coding': 14.285714285714286, 'Non-coding': 14.285714285714286,
        'Low complexity nucleotide': 14.285714285714286,
        'Read length was shorter than 3 * peptide k-mer size': 14.285714285714286,
        f'Low complexity peptide in {canonical_alphabet} alphabet': 14.285714285714286}
    true_percentages = {
        'All translations shorter than peptide k-mer size + 1': 1,
        'All translation frames have stop codons': 1, 'Coding': 1,
        'Non-coding': 1, 'Low complexity nucleotide': 1,
        'Read length was shorter than 3 * peptide k-mer size': 1,
        f'Low complexity peptide in {canonical_alphabet} alphabet': 1}
    assert test_counts == true_counts
    assert test_percentages == true_percentages


def test_generate_coding_summary(
        reads, data_folder,
        single_alphabet_ksize_true_scores):
    assemble_ss = AssembleSaveSummary(
        reads, True, True,
        'bloom_filter.nodegraph',
        "protein", 7, 0.5)
    test_summary = assemble_ss.generate_coding_summary(
        single_alphabet_ksize_true_scores)

    true_summary = {
        'input_files': [
            'SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled_n22.fq'],
        'jaccard_info': {'count': 17.0, 'mean': 0.2186899269511726,
                         'std': 0.37622482326071616, 'min': 0.0,
                         '25%': 0.0, '50%': 0.0625, '75%': 0.125,
                         'max': 1.0},
        'classification_value_counts': {
            'All translations shorter than peptide k-mer size + 1': 0,
            'All translation frames have stop codons': 3,
            'Coding': 3, 'Non-coding': 14,
            'Low complexity nucleotide': 0,
            'Read length was shorter than 3 * peptide k-mer size': 2,
            'Low complexity peptide in protein20 alphabet': 1},
        'classification_percentages': {
            'All translations shorter than peptide k-mer size + 1': 0.0,
            'All translation frames have stop codons': 13.043478260869565,
            'Coding': 13.043478260869565,
            'Non-coding': 60.869565217391305,
            'Low complexity nucleotide': 0.0,
            'Read length was shorter than 3 * peptide k-mer size': 8.695652173913043,
            'Low complexity peptide in protein20 alphabet': 4.3478260869565215},
        'histogram_n_coding_frames_per_read': {
            'Number of reads with 1 putative protein-coding translations': 3},
        'histogram_n_coding_frames_per_read_percentages': {
            'Number of reads with 1 putative protein-coding translations': 100.0},
        'peptide_bloom_filter': 'bloom_filter.nodegraph',
        'peptide_alphabet': 'protein', 'peptide_ksize': 7,
        'jaccard_threshold': 0.5}
    assert test_summary == true_summary


def test_maybe_write_csv(
        reads, single_alphabet_ksize_true_scores, true_scores_path):
    assemble_ss = AssembleSaveSummary(
        reads, true_scores_path, True,
        'bloom_filter.nodegraph',
        "protein", 7, 0.5)
    assemble_ss.maybe_write_csv(
        single_alphabet_ksize_true_scores)


def test_make_empty_coding_categories():
    assemble_ss = AssembleSaveSummary(
        ['nonexistent.fa'], True, True,
        'bloom_filter.nodegraph',
        "protein", 7, 0.5)
    test_coding_categories = {
        'All translations shorter than peptide k-mer size + 1': 0,
        'All translation frames have stop codons': 0,
        'Coding': 0,
        'Non-coding': 0,
        'Low complexity nucleotide': 0,
        'Read length was shorter than 3 * peptide k-mer size': 0,
        'Low complexity peptide in protein20 alphabet': 0}
    true_coding_categories = assemble_ss.make_empty_coding_categories()
    assert true_coding_categories == test_coding_categories
