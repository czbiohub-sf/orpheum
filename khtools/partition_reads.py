"""
partition_reads.py

Partition reads into coding, noncoding, and low-complexity bins
"""
import gzip
from io import StringIO
from pprint import pprint
import warnings

from Bio.Seq import Seq
from Bio import SeqIO
import click
from khmer import Nodegraph
import matplotlib.pyplot as plt
import screed
from sourmash._minhash import hash_murmur
from khmer.khmer_args import calculate_graphsize
import pandas as pd
from sourmash.logging import notify
import seaborn as sns
from sklearn import metrics
from khtools.compare_kmer_content import kmerize, hpize, dayhoffize
from tqdm import tqdm


# Import modified 'os' module with LC_LANG set so click doesn't complain.
# The '# noqa: F401' line prevents the linter from complaining about the unused
# import.
from .os_utils import os    # noqa: F401


DEFAULT_K = 32
DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = 1e10
DEFAULT_N_THREADS = 1
VALID_MOLECULES = 'protein', 'peptide', 'dayhoff', 'hydrophobic-polar'
DEFAULT_SEED = 42

def make_peptide_bloom_filter(peptide_fasta, peptide_ksize, n_tables=4,
                              molecule='protein',
                              tablesize=DEFAULT_MAX_TABLESIZE,
                              seed=DEFAULT_SEED):
    """Create a bloom filter out of peptide seuqences"""
    try:
        assert molecule in VALID_MOLECULES
    except AssertionError:
        raise ValueError(f"{molecule} is not a valid amino acid molecule " \
                          "type. Only {','.join(VALID_MOLECULES)} are " \
                          "accepted")

    peptide_graph = Nodegraph(peptide_ksize, tablesize, n_tables=n_tables,
                              seed=seed)

    for record in screed.open(peptide_fasta):
        if '*' in record['sequence']:
            continue
        sequence = record['sequence']
        if molecule == 'dayhoff':
            sequence = dayhoffize(sequence)
        elif molecule == 'hydrophobic-polar':
            sequence = hpize(sequence)
        kmers = kmerize(sequence, peptide_ksize)
        for kmer in kmers:
            # Convert the k-mer into an integer
            hashed = hash_murmur(kmer)

            # .add can take the hashed integer so we can hash the peptide
            #  kmer and add it directly
            peptide_graph.add(hashed)
    return peptide_graph


def three_frame_translation(seq, cds=False, debug=False):
    if debug:
        warning_filter = 'default'
    else:
        warning_filter = 'ignore'

    with warnings.catch_warnings():
        warnings.simplefilter(warning_filter)
        for frame in range(3):
            translation = seq[frame:].translate()
            yield translation


def three_frame_translation_no_stops(seq, cds=False, debug=False):
    return [t for t in three_frame_translation(seq, cds, debug)
            if '*' not in t]


def six_frame_translation_no_stops(seq, cds=False, debug=False):
    forward_translations = three_frame_translation_no_stops(seq, cds, debug)
    reverse_translations = three_frame_translation_no_stops(
        seq.reverse_complement(), cds, debug)
    return forward_translations + reverse_translations


@click.command()
@click.argument('reads')
@click.argument('peptides')
@click.option('--peptide-ksize', default=7,
                help="K-mer size of the peptide sequence to use")
@click.option("--long-reads", is_flag=True,
              help="If set, then only considers reading frames starting with "
                   "start codon (ATG) and ending in a stop codon "
                   "(TAG, TAA, TGA)")
@click.option("--verbose", is_flag=True,
              help="Print more output")
@click.option("--debug", is_flag=True,
                  help="Print developer debugger output, including warnings")
def cli(reads, peptides, peptide_ksize, long_reads=False, verbose=False):
    """

    Parameters
    ----------
    reads : str
        Sequence file of reads to filter
    peptides : str
        Sequence file of peptides
    peptide_ksize : int

    long_reads
    verbose

    Returns
    -------

    """
    peptide_graph = make_peptide_bloom_filter(peptides, peptide_ksize)

    for x in tqdm(range(count)):
        # note that colorama.init() doesn't need to be called for the colors
        # to work
        click.echo(click.style('Hello %s!' % name, fg=random.choice(COLORS)))
