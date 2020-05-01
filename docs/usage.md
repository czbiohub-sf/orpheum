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




## Best K-mer size for each alphabet

Emprically, a DNA k-mer of 21 and a protein k-mer of 7 is sufficient to distinguish transcriptomes. Using this "ideal" k-size of 7 for proteins (which have an alphabet size of 7), then if <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;|\Sigma|" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;|\Sigma|" title="|\Sigma|" /></a> represents the size of the k-mer alphabet, the ideal k-mer size for each alphabet size can be derived this way:

<a href="https://www.codecogs.com/eqnedit.php?latex=20^7&space;=&space;|\Sigma|^x\\&space;7&space;\log&space;20&space;=&space;x&space;\log&space;|\Sigma|\\&space;x&space;=&space;\frac{7&space;\log&space;20}{\log&space;|\Sigma|}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?20^7&space;=&space;|\Sigma|^x\\&space;7&space;\log&space;20&space;=&space;x&space;\log&space;|\Sigma|\\&space;x&space;=&space;\frac{7&space;\log&space;20}{\log&space;|\Sigma|}" title="20^7 = |\Sigma|^x\\ 7 \log 20 = x \log |\Sigma|\\ x = \frac{7 \log 20}{\log |\Sigma|}" /></a>

Below is a table of the raw float values for this, and if you take the ceiling (always round up) of each decimal into an integer:

| ksize_float        | alphabet_size | ksize | 
|--------------------|---------------|-------| 
| 30.253496664211536 | 2             | 31    | 
| 19.087831195025892 | 3             | 20    | 
| 15.126748332105768 | 4             | 16    | 
| 13.029471813027504 | 5             | 14    | 
| 11.703650113211072 | 6             | 12    | 
| 10.776512946940723 | 7             | 11    | 
| 10.084498888070513 | 8             | 11    | 
| 9.543915597512946  | 9             | 10    | 
| 9.107209969647867  | 10            | 10    | 
| 8.745221758749107  | 11            | 9     | 
| 8.438999475761797  | 12            | 9     | 
| 8.175649103509599  | 13            | 9     | 
| 7.946066832104447  | 14            | 8     | 
| 7.7436252497619265 | 15            | 8     | 
| 7.563374166052884  | 16            | 8     | 
| 7.401534359871294  | 17            | 8     | 
| 7.255165657355305  | 18            | 8     | 
| 7.1219427752632525 | 19            | 8     | 
| 7.0                | 20            | 7     | 
