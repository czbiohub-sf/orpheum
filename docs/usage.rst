========
Usage
========

To use Kmer-hashing tools in a project::

    import khtools

To create a bloom filter of sequences::

    khtools bloom-filter --molecule protein --peptide-ksize 7 --save-as Homo_sapiens.GRCh38.pep.subset.molecule-protein_ksize-7.bloomfilter.nodegraph Homo_sapiens.GRCh38.pep.subset.fa.gz

To partition reads into coding/noncoding bins using the bloom filter::

    khtools partition -- SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled.fq.gz Homo_sapiens.GRCh38.pep.all.fa.gz

To create the bloom filter and partition the reads in one step::

    khtools partition  ~/code/kmer-hashing/extract_kmers/test-data/SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled.fq.gz ~/Downloads/Homo_sapiens.GRCh38.pep.all.fa.gz
