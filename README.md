Kmer-hashing tools
================================

[![image](https://img.shields.io/travis/khtools.svg)](https://travis-ci.org/khtools)


[![codecov](https://codecov.io/gh/khtools/branch/master/graph/badge.svg)](https://codecov.io/gh/khtools)

[![image](https://img.shields.io/pypi/v/%7B%7B%20cookiecutter.repo_name%20%7D%7D.svg)](https://pypi.python.org/pypi/%7B%7B%20cookiecutter.repo_name%20%7D%7D)


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

The "coding score" is obtained by the
[Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index) between any of the
 six translated frames of the read and the peptide database.

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

