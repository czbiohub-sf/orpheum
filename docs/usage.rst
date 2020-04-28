========
Usage
========

To use sencha in a project::

    import sencha

To create a protein index::

    sencha index --molecule protein --peptide-ksize 7 --save-as Homo_sapiens.GRCh38.pep.subset.molecule-protein_ksize-7.bloomfilter.nodegraph Homo_sapiens.GRCh38.pep.subset.fa.gz

To translate RNA-seq reads into coding peptides using the protein index::

    sencha translate -- SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled.fq.gz Homo_sapiens.GRCh38.pep.all.fa.gz

To create the index and translate the reads in one step::

    sencha partition  ~/code/kmer-hashing/extract_kmers/test-data/SRR306838_GSM752691_hsa_br_F_1_trimmed_subsampled.fq.gz ~/Downloads/Homo_sapiens.GRCh38.pep.all.fa.gz
