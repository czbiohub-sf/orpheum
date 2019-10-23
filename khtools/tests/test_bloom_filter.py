import math


def test_make_peptide_bloom_filter(peptide_fasta, molecule, peptide_ksize):
    from khtools.bloom_filter import make_peptide_bloom_filter

    test = make_peptide_bloom_filter(peptide_fasta, peptide_ksize, molecule,
                                     tablesize=1e4)
    TRUE_N_UNIQUE_KMERS = {("protein", 7): 20717,
                           ("dayhoff", 7): 20669,
                           ("dayhoff", 12): 20712,
                           ("hydrophobic-polar", 21): 20712,
                           ("hydrophobic-polar", 7): 170}
    true_n_unique_kmers = TRUE_N_UNIQUE_KMERS[(molecule, peptide_ksize)]

    # For now, assert that the number of kmers is within 5% of the true value
    assert test.n_unique_kmers() > true_n_unique_kmers * 0.95
    assert test.n_unique_kmers() < true_n_unique_kmers * 1.05

