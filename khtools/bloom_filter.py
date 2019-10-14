import os

import click
from khmer import Nodegraph
import screed
from sourmash._minhash import hash_murmur

from khtools.compare_kmer_content import kmerize
from khtools.sequence_encodings import encode_peptide


DEFAULT_MAX_TABLESIZE = 1e10


def make_peptide_bloom_filter(peptide_fasta, peptide_ksize, molecule='protein',
                              n_tables=4, tablesize=DEFAULT_MAX_TABLESIZE):
    """Create a bloom filter out of peptide sequences"""
    peptide_graph = Nodegraph(peptide_ksize, tablesize, n_tables=n_tables)

    with screed.open(peptide_fasta) as records:
        for record in records:
            if '*' in record['sequence']:
                continue
            sequence = encode_peptide(record['sequence'], molecule)
            kmers = kmerize(sequence, peptide_ksize)
            for kmer in kmers:
                # Convert the k-mer into an integer
                hashed = hash_murmur(kmer)

                # .add can take the hashed integer so we can hash the peptide
                #  kmer and add it directly
                peptide_graph.add(hashed)
    return peptide_graph


def maybe_make_peptide_bloom_filter(peptides, peptide_ksize,
                                    molecule,
                                    peptides_are_bloom_filter):
    if peptides_are_bloom_filter:
        print(f"Loading existing bloom filter from {peptides} and making " \
               "sure the ksizes match")
        peptide_graph = Nodegraph.load(peptides)
        assert peptide_ksize == peptide_graph.ksize
    else:
        print(f"Creating peptide bloom filter with file: {peptides} using " \
               f"ksize: {peptide_ksize} and molecule: {molecule} ...")
        peptide_graph = make_peptide_bloom_filter(peptides, peptide_ksize,
                                                  molecule=molecule)
    return peptide_graph



def maybe_save_peptide_bloom_filter(peptides, peptide_graph,
                                    molecule, ksize,
                                    save_peptide_bloom_filter):
    if save_peptide_bloom_filter:

        if isinstance(save_peptide_bloom_filter, str):
            filename = save_peptide_bloom_filter
            peptide_graph.save(save_peptide_bloom_filter)
        else:
            suffix = f'.molecule-{molecule}_ksize-{ksize}.bloomfilter.nodegraph'
            filename = os.path.splitext(peptides)[0] + suffix

        print(f"Writing peptide bloom filter to {filename}")
        peptide_graph.save(filename)
        print("\tDone!")


@click.command()
@click.argument('peptides')
@click.option('--peptide-ksize', default=7,
                help="K-mer size of the peptide sequence to use. Default: 7")
@click.option('--molecule', default='protein',
              help="The type of amino acid encoding to use. Default is "
                   "'protein', but 'dayhoff' or 'hydrophobic-polar' can be "
                   "used")
@click.option('--save-as', default=None,
              help='If provided, save peptide bloom filter as this filename. '
                   'Otherwise, add ksize and molecule name to input filename.')
def cli(peptides, peptide_ksize=7, molecule='protein', save_as=None):
    """

    \b
    Parameters
    ----------
    reads : str
        Sequence file of reads to filter
    peptides : str
        Sequence file of peptides
    peptide_ksize : int
        Number of characters in amino acid words
    long_reads
    verbose

    \b
    Returns
    -------

    """
    # \b above prevents rewrapping of paragraph
    peptide_graph = make_peptide_bloom_filter(peptides, peptide_ksize,
                                              molecule)
    click.echo("\tDone!")

    save_peptide_bloom_filter = save_as if save_as is not None else True
    maybe_save_peptide_bloom_filter(peptides, peptide_graph,
                                    molecule, peptide_ksize,
                                    save_peptide_bloom_filter=save_peptide_bloom_filter)

