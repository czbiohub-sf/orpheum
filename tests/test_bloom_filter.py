

def test_make_peptide_bloom_filter(peptide_fasta, molecule, peptide_ksize):
    from khtools.bloom_filter import make_peptide_bloom_filter

    test = make_peptide_bloom_filter(peptide_fasta,
                                     peptide_ksize,
                                     molecule,
                                     tablesize=1e6)
    TRUE_N_UNIQUE_KMERS = {
        ("protein", 7): 506352,
        ("dayhoff", 7): 99863,
        ("dayhoff", 11): 472197,
        ("dayhoff", 12): 488469,
        ("hydrophobic-polar", 21): 434810,
        ("hydrophobic-polar", 7): 170
    }
    true_n_unique_kmers = TRUE_N_UNIQUE_KMERS[(molecule, peptide_ksize)]

    # For now, assert that the number of kmers is within 5% of the true value
    assert test.n_unique_kmers() > true_n_unique_kmers * 0.95
    assert test.n_unique_kmers() < true_n_unique_kmers * 1.05
