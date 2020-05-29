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

#### Save the "coding scores" to a csv

The "coding score" of each read is calculated by translating each read in six
frames, then is calculatating the
[Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index) between any of the
six translated frames of the read and the peptide database. The final coding
score is the maximum Jaccard index across all reading frames. If you'd like to
see the coding scores for all reads, use the `--csv` flag.

```
sencha translate --csv coding_scores.csv reference-proteome.fa.gz *.fastq.gz > coding_peptides.fasta
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

## Checks and balances to ensure a worry-free translation

Some of the concerns with translating reads into protein is what is the false positive rate of the translation? How can I be sure that only "real" reading frames are being passed through and not just all possible k-mers?

There are two main things contributing to the false positive rate:

1. The k-mer size and alphabet
2. The false positive of the bloom filter

### How big of a k-mer size should I use for my alphabet?

To find true protein-coding translations, `sencha` relies on having a large "negative space" of unobserved k-mers which represent not real protein coding k-mers. Empirically, the ratio of observed k-mers to theoretical k-mers should be about 1e-4 to ensure that there are enough unobserved k-mers to consider as non-coding.

The default k-mer sizes for the various alphabets and their alphabet sizes (sigma):

- Protein (sigma = 20): 9
- Dayhoff (sigma = 6): 15
- Hydrophobic-Polar (sigma = 2): 39

These numbers came from empirical testing using the UniProt/SwissProt manually reviewed database of mammalian protein sequences on primate data, and a protein k-size of 9 created "reasonable" results. Meaning, there were few reads that had multiple translations per read, and especially not 3+ translations per read. Extrapolating to other k-mer sizes, one can use log math to solve for `k` in this equation: `20^9 = sigma^k`, where `sigma` is the size of the alphabet.

Now, `sencha` will throw an error if it senses that the size k-mers from the provided reference proteome is too large relative to the total theoretical k-mers. Here is what that error looks like:

```
 ✘  Fri 29 May - 14:28  ~/code/sencha   origin ☊ olgabot/index-updates ↑1 2☀ 2● 
  sencha index --molecule protein --peptide-ksize 7 tests/data/index/Homo_sapiens.GRCh38.pep.subset.fa.gz
2662it [00:05, 510.07it/s]
Traceback (most recent call last):
  File "/anaconda3/envs/sencha-master/bin/sencha", line 11, in <module>
    load_entry_point('sencha', 'console_scripts', 'sencha')()
  File "/anaconda3/envs/sencha-master/lib/python3.7/site-packages/click/core.py", line 829, in __call__
    return self.main(*args, **kwargs)
  File "/anaconda3/envs/sencha-master/lib/python3.7/site-packages/click/core.py", line 782, in main
    rv = self.invoke(ctx)
  File "/anaconda3/envs/sencha-master/lib/python3.7/site-packages/click/core.py", line 1259, in invoke
    return _process_result(sub_ctx.command.invoke(sub_ctx))
  File "/anaconda3/envs/sencha-master/lib/python3.7/site-packages/click/core.py", line 1066, in invoke
    return ctx.invoke(self.callback, **ctx.params)
  File "/anaconda3/envs/sencha-master/lib/python3.7/site-packages/click/core.py", line 610, in invoke
    return callback(*args, **kwargs)
  File "/Users/olgabot/code/sencha/sencha/index.py", line 322, in cli
    max_observed_fraction=max_observed_fraction,
  File "/Users/olgabot/code/sencha/sencha/index.py", line 121, in make_peptide_index
    check_kmer_occupancy(max_observed_fraction, molecule, peptide_index, peptide_ksize)
  File "/Users/olgabot/code/sencha/sencha/index.py", line 133, in check_kmer_occupancy
    f"The number of observed length {peptide_ksize} k-mers compared to the "
ValueError: The number of observed length 7 k-mers compared to the possible theoretical k-mers is 509247 / 1280000000 = 3.98e-04 which is greater than the maximum observed fraction threshold, 1.00e-04. This doesn't leave enough 'negative space' of non-observed protein k-mers for room for predicting true protein-coding sequence, which relies on seeing which protein k-mers are *not* present in the observed data. Please increase the k-mer size.
```

### How big of a table size should I use?

The default "tablesize" is 1e8, which has been sufficient for UniProt/SwissProt reference proteome from Mammalia. However, if you are using a more complex dataset with many more k-mers, you may want to increase the tablesize to avoid false positive collisions in the bloom filter index.

With a `--tablesize` that is too small, the error looks like this:

```
(sencha-master) 
 ✘  Fri 29 May - 14:29  ~/code/sencha   origin ☊ olgabot/index-updates ↑1 2☀ 2● 
  sencha index --tablesize 1e4 tests/data/index/Homo_sapiens.GRCh38.pep.subset.fa.gz            
2662it [00:03, 667.66it/s]
**
** ERROR: the graph structure is too small for 
** this data set.  Increase data structure size
** with --max_memory_usage/-M.
**
** Do not use these results!!
**
** (estimated false positive rate of 1.013; max recommended 0.200)
**
Traceback (most recent call last):
  File "/anaconda3/envs/sencha-master/bin/sencha", line 11, in <module>
    load_entry_point('sencha', 'console_scripts', 'sencha')()
  File "/anaconda3/envs/sencha-master/lib/python3.7/site-packages/click/core.py", line 829, in __call__
    return self.main(*args, **kwargs)
  File "/anaconda3/envs/sencha-master/lib/python3.7/site-packages/click/core.py", line 782, in main
    rv = self.invoke(ctx)
  File "/anaconda3/envs/sencha-master/lib/python3.7/site-packages/click/core.py", line 1259, in invoke
    return _process_result(sub_ctx.command.invoke(sub_ctx))
  File "/anaconda3/envs/sencha-master/lib/python3.7/site-packages/click/core.py", line 1066, in invoke
    return ctx.invoke(self.callback, **ctx.params)
  File "/anaconda3/envs/sencha-master/lib/python3.7/site-packages/click/core.py", line 610, in invoke
    return callback(*args, **kwargs)
  File "/Users/olgabot/code/sencha/sencha/index.py", line 322, in cli
    max_observed_fraction=max_observed_fraction,
  File "/Users/olgabot/code/sencha/sencha/index.py", line 119, in make_peptide_index
    check_collisions(peptide_index, tablesize)
  File "/Users/olgabot/code/sencha/sencha/index.py", line 150, in check_collisions
    f"The false positive rate in the bloom filter index is "
ValueError: The false positive rate in the bloom filter index is 1.0129382731984724, which is greater than than the recommended maximum of 0.2. The current table size is 1.0e+04, please increase by an order of magnitude and rerun.

```