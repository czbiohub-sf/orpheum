import os

import click
import khmer
import screed
from sourmash._minhash import hash_murmur
from tqdm import tqdm

from khtools.compare_kmer_content import kmerize
from khtools.sequence_encodings import encode_peptide, VALID_PEPTIDE_MOLECULES

# khmer Nodegraph features
DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = int(1e8)

# Default k-mer sizes for different alphabets
DEFAULT_PROTEIN_KSIZE = 7
DEFAULT_DAYHOFF_KSIZE = 11
DEFAULT_HP_KSIZE = 21


def load_nodegraph(*args, **kwargs):
    try:
        # khmer 2.1.1
        return khmer.load_nodegraph(*args, **kwargs)
    except AttributeError:
        # khmer 3+/master branch
        return khmer.Nodegraph.load(*args, **kwargs)


# Cribbed from https://click.palletsprojects.com/en/7.x/parameters/
class BasedIntParamType(click.ParamType):
    name = "integer"

    def convert(self, value, param, ctx):
        try:
            if isinstance(value, int):
                return value
            if 'e' in value:
                sigfig, exponent = value.split('e')
                sigfig = float(sigfig)
                exponent = int(exponent)
                return int(sigfig * 10 ** exponent)
            return int(value, 10)
        except TypeError:
            self.fail(
                "expected string for int() conversion, got "
                f"{value!r} of type {type(value).__name__}",
                param,
                ctx,
            )
        except ValueError:
            self.fail(f"{value!r} is not a valid integer", param, ctx)


BASED_INT = BasedIntParamType()


def make_peptide_bloom_filter(peptide_fasta,
                              peptide_ksize,
                              molecule,
                              n_tables=DEFAULT_N_TABLES,
                              tablesize=DEFAULT_MAX_TABLESIZE):
    """Create a bloom filter out of peptide sequences"""
    peptide_bloom_filter = khmer.Nodegraph(peptide_ksize,
                                           tablesize,
                                           n_tables=n_tables)

    with screed.open(peptide_fasta) as records:
        for record in tqdm(records):
            if '*' in record['sequence']:
                continue
            sequence = encode_peptide(record['sequence'], molecule)
            try:
                kmers = kmerize(sequence, peptide_ksize)
                for kmer in kmers:
                    # Convert the k-mer into an integer
                    hashed = hash_murmur(kmer)

                    # .add can take the hashed integer so we can hash the
                    #  peptide kmer and add it directly
                    peptide_bloom_filter.add(hashed)
            except ValueError:
                # Sequence length is smaller than k-mer size
                continue
    return peptide_bloom_filter


def maybe_make_peptide_bloom_filter(peptides, peptide_ksize, molecule,
                                    peptides_are_bloom_filter,
                                    n_tables=DEFAULT_N_TABLES,
                                    tablesize=DEFAULT_MAX_TABLESIZE):
    if peptides_are_bloom_filter:
        click.echo(
            f"Loading existing bloom filter from {peptides} and "
            f"making sure the ksizes match",
            err=True)
        peptide_bloom_filter = load_nodegraph(peptides)
        if peptide_ksize is not None:
            try:
                assert peptide_ksize == peptide_bloom_filter.ksize()
            except AssertionError:
                raise ValueError(f"Given peptide ksize ({peptide_ksize}) and "
                                 f"ksize found in bloom filter "
                                 f"({peptide_bloom_filter.ksize()}) are not"
                                 f"equal")
    else:
        click.echo(
            f"Creating peptide bloom filter with file: {peptides}\n"
            f"Using ksize: {peptide_ksize} and molecule: {molecule} "
            f"...",
            err=True)
        peptide_bloom_filter = make_peptide_bloom_filter(
            peptides, peptide_ksize, molecule=molecule,
            n_tables=n_tables, tablesize=tablesize)
    return peptide_bloom_filter


def maybe_save_peptide_bloom_filter(peptides, peptide_bloom_filter, molecule,
                                    save_peptide_bloom_filter):
    if save_peptide_bloom_filter:
        ksize = peptide_bloom_filter.ksize()

        if isinstance(save_peptide_bloom_filter, str):
            filename = save_peptide_bloom_filter
            peptide_bloom_filter.save(save_peptide_bloom_filter)
        else:
            suffix = f'.molecule-{molecule}_ksize-{ksize}.bloomfilter.' \
                     f'nodegraph'
            filename = os.path.splitext(peptides)[0] + suffix

        click.echo(f"Writing peptide bloom filter to {filename}", err=True)
        peptide_bloom_filter.save(filename)
        click.echo("\tDone!", err=True)


@click.command()
@click.argument('peptides')
@click.option('--peptide-ksize',
              default=None, type=int,
              help="K-mer size of the peptide sequence to use. Defaults for"
              " different molecules are, "
              f"protein: {DEFAULT_PROTEIN_KSIZE}"
              f", dayhoff: {DEFAULT_DAYHOFF_KSIZE},"
              f" hydrophobic-polar: {DEFAULT_HP_KSIZE}")
@click.option('--molecule',
              default='protein',
              help="The type of amino acid encoding to use. Default is "
              "'protein', but 'dayhoff' or 'hydrophobic-polar' can be "
              "used")
@click.option('--save-as',
              default=None,
              help='If provided, save peptide bloom filter as this filename. '
              'Otherwise, add ksize and molecule name to input filename.')
@click.option('--tablesize', type=BASED_INT,
              default="1e8",
              help='Size of the bloom filter table to use')
@click.option('--n-tables', type=int,
              default=DEFAULT_N_TABLES,
              help='Size of the bloom filter table to use')
def cli(peptides, peptide_ksize=None, molecule='protein', save_as=None,
        tablesize=DEFAULT_MAX_TABLESIZE, n_tables=DEFAULT_N_TABLES):
    """Make a peptide bloom filter for your peptides

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
    peptide_ksize = get_peptide_ksize(molecule, peptide_ksize)
    peptide_bloom_filter = make_peptide_bloom_filter(peptides, peptide_ksize,
                                                     molecule,
                                                     n_tables=n_tables,
                                                     tablesize=tablesize)
    click.echo("\tDone!", err=True)

    save_peptide_bloom_filter = save_as if save_as is not None else True
    maybe_save_peptide_bloom_filter(
        peptides,
        peptide_bloom_filter,
        molecule,
        save_peptide_bloom_filter=save_peptide_bloom_filter)


def get_peptide_ksize(molecule, peptide_ksize):
    if molecule not in VALID_PEPTIDE_MOLECULES:
        raise ValueError(f"{molecule} is not a valid protein encoding! "
                         f"Only one of 'protein', 'hydrophobic-polar', or"
                         f" 'dayhoff' can be specified")

    if peptide_ksize is None:
        if molecule == 'protein':
            peptide_ksize = DEFAULT_PROTEIN_KSIZE
        elif molecule == 'dayhoff':
            peptide_ksize = DEFAULT_DAYHOFF_KSIZE
        elif molecule == 'hydrophobic-polar' or molecule == 'hp':
            peptide_ksize = DEFAULT_HP_KSIZE
    return peptide_ksize
