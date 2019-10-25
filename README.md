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
 score is the maximum Jaccard index across all reading frames.

```
khtools extract_coding --csv coding_scores.csv peptides.fa.gz *.fastq.gz > coding_peptides.fasta
```


#### Save the coding nucleotides to a fasta

```
khtools extract_coding --coding-nucleotide-fasta coding_nucleotides.fasta peptides.fa.gz *.fastq.gz > coding_peptides.fasta
```

#### Save the *non*-coding nucleotides to a fasta

```
khtools extract_coding --noncoding-nucleotide-fasta noncoding_nucleotides.fasta peptides.fa.gz *.fastq.gz > coding_peptides.fasta
```


Features
--------

-   TODO

