Kmer-hashing tools
================================

[![image](https://img.shields.io/travis/czbiohub/kh-tools.svg)](https://travis-ci.com/czbiohub/kh-tools)
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
would be considered low complexity.

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
or `CCCCCCCCCCCCCCCCCCCCCCCCCCCCC` (3'5' Frame 2).

```
khtools extract_coding --low-complexity-peptides-fasta low_complexity_peptides.fasta peptides.fa.gz *.fastq.gz > coding_peptides.fasta
```

