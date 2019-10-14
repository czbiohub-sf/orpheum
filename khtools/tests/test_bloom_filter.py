


def test_make_peptide_bloom_filter(peptide_fasta, molecule, peptide_ksize):
    from khtools.bloom_filter import make_peptide_bloom_filter

    test = make_peptide_bloom_filter(peptide_fasta, peptide_ksize, molecule,
                                     tablesize=1e6)
    TRUE_N_UNIQUE_KMERS = {("protein", 7): 506352,
                           ("protein", 8): 512787,
                           ("dayhoff", 7): 99863,
                           ("dayhoff", 8): 223991,
                           ("hydrophobic-polar", 7): 170,
                           ("hydrophobic-polar", 8): 317}
    true_n_unique_kmers = TRUE_N_UNIQUE_KMERS[(molecule, peptide_ksize)]
    assert test.n_unique_kmers() == true_n_unique_kmers
    filename = f'/Users/olgabot/code/kmer-hashing/kh-tools/khtools/tests/data/bloom_filter/Homo_sapiens.GRCh38.pep.subset.molecule-{protein}_ksize-{peptide_ksize}.bloomfilter.nodegraph'
    test.save(filename)

