sencha
================================
![Tests](https://github.com/czbiohub/sencha/workflows/Pytest/badge.svg)
![Linting](https://github.com/czbiohub/sencha/workflows/Lint%20with%20flake8/badge.svg)
[![codecov](https://codecov.io/gh/czbiohub/sencha/branch/master/graph/badge.svg)](https://codecov.io/gh/czbiohub/sencha)

What is sencha?
-------------------------------------

Sencha is a Python package for directly translating RNA-seq reads into coding protein sequence.

-   Free software: MIT license
-   Documentation: https://czbiohub.github.io/sencha

The name is inspired from the naming pattern of [sourmash](https://github.com/dib-lab/sourmash)
combined with [@olgabot](https://github.com/olgabot/)'s love of tea.
([Sencha](https://en.wikipedia.org/wiki/Sencha) is a Japanese green tea.)

Installation
------------

The package can be installed from PyPI using `pip` here:

```
pip install sencha
```

### Developmental install

To install this code and play around with the code locally, clone this github repository and use `pip` to install:

```
git clone https://github.com/czbiohub/sencha.git
cd sencha

# The "." means "install *this*, the folder where I am now"
pip install .
```

Usage
-----

### Extract likely protein-coding reads from sequencing data

A reference proteome *must* be supplied as the first argument.

```
sencha translate reference-proteome.fa.gz *.fastq.gz > coding_peptides.fasta
```

#### Save the "coding scores" to a csv or parquet file

The "coding score" of each read is calculated by translating each read in six
frames, then is calculatating the
[Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index) between any of the
six translated frames of the read and the peptide database. The final coding
score is the maximum Jaccard index across all reading frames. If you'd like to
see the coding scores for all reads, use the `--csv` flag or `--parquet` flag.

csv:
```
sencha translate --csv coding_scores.csv reference-proteome.fa.gz *.fastq.gz > coding_peptides.fasta
```

parquet:
```
sencha translate --parquet coding_scores.parquet reference-proteome.fa.gz *.fastq.gz > coding_peptides.fasta
```

#### Save the coding nucleotides to a fasta

By default, only the coding *peptides* are output. If you'd like to also output
the underlying *nucleotide* sequence, then use the flag `--coding-nucleotide-fasta`

```
sencha translate --coding-nucleotide-fasta coding_nucleotides.fasta reference-proteome.fa.gz *.fastq.gz > coding_peptides.fasta
```

#### Save the *non*-coding nucleotides to a fasta

To see the sequence of reads which were deemed non-coding, use the flag
`--noncoding-nucleotide-fasta`.

```
sencha translate --noncoding-nucleotide-fasta noncoding_nucleotides.fasta reference-proteome.fa.gz *.fastq.gz > coding_peptides.fasta
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
sencha translate --low-complexity-nucleotide-fasta low_complexity_nucleotides.fasta reference-proteome.fa.gz *.fastq.gz > coding_peptides.fasta
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
sencha translate --low-complexity-peptides-fasta low_complexity_peptides.fasta reference-proteome.fa.gz *.fastq.gz > coding_peptides.fasta
```

