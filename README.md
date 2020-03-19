Kmer-hashing tools
================================
![Tests](https://github.com/czbiohub/kh-tools/workflows/Pytest/badge.svg)
![Linting](https://github.com/czbiohub/kh-tools/workflows/Lint%20with%20flake8/badge.svg)
[![codecov](https://codecov.io/gh/czbiohub/kh-tools/branch/master/graph/badge.svg)](https://codecov.io/gh/czbiohub/kh-tools)

What is khtools?
-------------------------------------

Kmer hashing tools contains data cleaning and visualization code for analyzing kmer-hashing similarity matrices

-   Free software: MIT license
-   Documentation: <https://>czbiohub.github.io/khtools

Installation
------------

To install this code, clone this github repository and use pip to install

```
git clone <https://github.com/>czbiohub/khtools.git
cd khtools

# The "." means "install *this*, the folder where I am now"
pip install .
```

Usage
-----

### Extract likely protein-coding reads from sequencing data


```
khtools extract_coding peptides.fa.gz *.fastq.gz > coding_peptides.fasta
```

#### Save the "coding scores" to a csv

The "coding score" of each read is calculated by translating each read in six
frames, then is calculatating the
[Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index) between any of the
six translated frames of the read and the peptide database. The final coding
score is the maximum Jaccard index across all reading frames. If you'd like to
see the coding scores for all reads, use the `--csv` flag.

```
khtools extract_coding --csv coding_scores.csv peptides.fa.gz *.fastq.gz > coding_peptides.fasta
```


#### Save the coding nucleotides to a fasta

By default, only the coding *peptides* are output. If you'd like to also output
the underlying *nucleotide* sequence, then use the flag `--coding-nucleotide-fasta`

```
khtools extract_coding --coding-nucleotide-fasta coding_nucleotides.fasta peptides.fa.gz *.fastq.gz > coding_peptides.fasta
```

#### Save the *non*-coding nucleotides to a fasta

To see the sequence of reads which were deemed non-coding, use the flag
`--noncoding-nucleotide-fasta`.

```
khtools extract_coding --noncoding-nucleotide-fasta noncoding_nucleotides.fasta peptides.fa.gz *.fastq.gz > coding_peptides.fasta
```

#### Save the low complexity nucleotides to a fasta

To see the sequence of reads found to have too low complexity of nucleotide
sequence to evaluate, use the flag `--low-complexity-nucleotide-fasta`. Low
complexity is determined by the same method as the read trimmer
[fastp](https://github.com/OpenGene/fastp) in which we calculate what
percentage of the sequence has consecutive runs of the same base,
or mathematically, how often `seq[i] = seq[i+1]`. The default threshold is
`0.3`. As an example, the sequence `CCCCCCCCCACCACCACCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCACACACCCCCAACACCC`
would be considered low complexity. While this sequence has many nucleotide
k-mers, it is likely a result of a sequencing error and we ignore it.

```
khtools extract_coding --low-complexity-nucleotide-fasta low_complexity_nucleotides.fasta peptides.fa.gz *.fastq.gz > coding_peptides.fasta
```

#### Save the low complexity peptides to a fasta

Even if the nucleotide sequence may pass the complexity filter, the peptide
sequence may still be low complexity. As an example, all translated frames of
the sequence
`CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG`
would be considered low complexity, as it translates to either
`QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ` (5'3' Frame 1),
`SSSSSSSSSSSSSSSSSSSSSSSSSSSSS` (5'3' Frame 2),
`AAAAAAAAAAAAAAAAAAAAAAAAAAAAA` (5'3' Frame 3 and 3'5' Frame 3),
`LLLLLLLLLLLLLLLLLLLLLLLLLLLLLL` (3'5' Frame 1),
or `CCCCCCCCCCCCCCCCCCCCCCCCCCCCC` (3'5' Frame 2). As these sequences have few
k-mers and are difficult to assess for how "coding" they are, we ignore them.
Unlike for nucleotides where we look at runs of consecutive bases, we require
the translated peptide to contain greater than `(L - k + 1)/2` k-mers, where
`L` is the length of the sequence and `k` is the k-mer size. To save the
sequence of low-complexity peptides to a fasta, use the flag
`--low-complexity-peptides-fasta`.

```
khtools extract_coding --low-complexity-peptides-fasta low_complexity_peptides.fasta peptides.fa.gz *.fastq.gz > coding_peptides.fasta
```



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
