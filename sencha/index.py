import math
import os

import click
import khmer
import screed
from sourmash._minhash import hash_murmur
from tqdm import tqdm

from sencha.compare_kmer_content import kmerize
from sencha.sequence_encodings import encode_peptide, BEST_KSIZES, ALPHABET_SIZES
import sencha.constants_index as constants_index
from sencha.log_utils import get_logger

logger = get_logger(__file__)


def per_translation_false_positive_rate(
    n_kmers_in_translation,
    n_total_kmers=1e7,
    n_hash_functions=constants_index.DEFAULT_N_TABLES,
    tablesize=constants_index.DEFAULT_MAX_TABLESIZE,
):
    exponent = -n_hash_functions * n_total_kmers / tablesize
    print(f"exponent: {exponent}")

    # Probability that a single k-mer is randomly in the data
    # per_kmer_fpr = math.pow(1 - math.exp(exponent), n_hash_functions)

    # Use built-in `exp1m` = exp - 1
    # - (exp - 1) = 1 - exp
    per_kmer_fpr = math.pow(-math.expm1(exponent), n_hash_functions)
    print(f"per kmer false positive rate: {per_kmer_fpr}")

    # Probability that the number of k-mers in the read are all false positives
    per_read_fpr = math.pow(per_kmer_fpr, n_kmers_in_translation)
    return per_read_fpr


def per_read_false_positive_coding_rate(
    read_length,
    peptide_ksize,
    n_total_kmers=4e7,
    alphabet="protein",
    n_hash_functions=constants_index.DEFAULT_N_TABLES,
    tablesize=constants_index.DEFAULT_MAX_TABLESIZE,
):
    """Compute the false positive rate that a translated k-mer randomly
    appears in the database"""

    alphabet_size = ALPHABET_SIZES[alphabet]
    n_kmers_of_size = min(n_total_kmers, alphabet_size ** peptide_ksize)

    false_positive_rate = 0
    for frame in range(3):
        translated_length = math.floor((read_length - frame) / 3)
        n_kmers_in_translation = translated_length - peptide_ksize + 1
        frame_fpr = per_translation_false_positive_rate(
            n_kmers_in_translation,
            n_kmers_of_size,
            n_hash_functions=n_hash_functions,
            tablesize=tablesize,
        )
        # multiply by two for both forward and reverse translation frames
        frame_fpr = 2 * frame_fpr
        false_positive_rate += frame_fpr
    return false_positive_rate


def load_nodegraph(*args, **kwargs):
    """Wrapper to load khmer-style bloom filter called a 'nodegraph'"""
    try:
        # khmer 2.1.1
        return khmer.load_nodegraph(*args, **kwargs)
    except AttributeError:
        # khmer 3+/master branch
        return khmer.Nodegraph.load(*args, **kwargs)


def maybe_read_peptide_file(peptide_file):
    records = []
    try:
        records = screed.open(peptide_file)
    except ValueError:
        logger.info(
            f"File {peptide_file} doesn't seem to be a sequence file, skipping. \n"
            f"..."
        )
    return records


def make_peptide_bloom_filter(
    peptides,
    peptide_ksize,
    molecule,
    n_tables=constants_index.DEFAULT_N_TABLES,
    tablesize=constants_index.DEFAULT_MAX_TABLESIZE,
    index_dir=None,
):
    """Create a bloom filter out of peptide sequences"""
    peptide_bloom_filter = khmer.Nodegraph(peptide_ksize, tablesize, n_tables=n_tables)
    if not index_dir:
        # if not a directory, should be single file. Convert to list to use for loop below.
        peptides = [peptides]
    for peptide_fasta in peptides:
        records = maybe_read_peptide_file(peptide_fasta)
        for record in tqdm(records):
            if "*" in record["sequence"]:
                continue
            sequence = encode_peptide(record["sequence"], molecule)
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


def make_peptide_set(peptide_fasta_files, peptide_ksize, molecule):
    """Create a python set out of peptide sequence k-mers

    For comparing to the bloom filter in storage and performance
    """
    peptide_set = set([])
    for peptide_fasta in peptide_fasta_files:
        records = maybe_read_peptide_file(peptide_fasta)
        for record in tqdm(records):
            if "*" in record["sequence"]:
                continue
            sequence = encode_peptide(record["sequence"], molecule)
            try:
                kmers = kmerize(sequence, peptide_ksize)
                peptide_set.update(kmers)
            except ValueError:
                # Sequence length is smaller than k-mer size
                continue
    return peptide_set


def maybe_make_peptide_bloom_filter(
    peptides,
    peptide_ksize,
    molecule,
    peptides_are_bloom_filter,
    n_tables=constants_index.DEFAULT_N_TABLES,
    tablesize=constants_index.DEFAULT_MAX_TABLESIZE,
    index_dir=None,
):
    if peptides_are_bloom_filter:
        logger.info(
            f"Loading existing bloom filter from {peptides} and "
            f"making sure the ksizes match"
        )
        peptide_bloom_filter = load_nodegraph(peptides)
        if peptide_ksize is not None:
            try:
                assert peptide_ksize == peptide_bloom_filter.ksize()
            except AssertionError:
                raise ValueError(
                    f"Given peptide ksize ({peptide_ksize}) and "
                    f"ksize found in bloom filter "
                    f"({peptide_bloom_filter.ksize()}) are not"
                    f"equal"
                )
    else:
        peptide_ksize = get_peptide_ksize(molecule, peptide_ksize)
        if not index_dir:
            logger.info(
                f"Creating peptide bloom filter with file: {peptides}\n"
                f"Using ksize: {peptide_ksize} and alphabet: {molecule} "
                f"..."
            )
        else:
            logger.info(
                f"Creating peptide bloom filter from files in directory: {index_dir}\n"
                f"Using ksize: {peptide_ksize} and alphabet: {molecule} "
                f"..."
            )

        peptide_bloom_filter = make_peptide_bloom_filter(
            peptides,
            peptide_ksize,
            molecule=molecule,
            n_tables=n_tables,
            tablesize=tablesize,
            index_dir=index_dir,
        )
    return peptide_bloom_filter


def maybe_save_peptide_bloom_filter(
    peptides, peptide_bloom_filter, molecule, save_peptide_bloom_filter, index_dir=None
):
    if save_peptide_bloom_filter:
        ksize = peptide_bloom_filter.ksize()

        if isinstance(save_peptide_bloom_filter, str):
            filename = save_peptide_bloom_filter
            peptide_bloom_filter.save(save_peptide_bloom_filter)
        else:
            suffix = f".alphabet-{molecule}_ksize-{ksize}.bloomfilter.nodegraph"
            if not index_dir:
                filename = os.path.splitext(peptides)[0] + suffix
            else:
                basename = os.path.basename(
                    index_dir
                )  # user index dir name as basename
                filename = os.path.join(index_dir, basename + suffix)

        logger.info(f"Writing peptide bloom filter to {filename}")
        peptide_bloom_filter.save(filename)
        logger.info("\tDone!")
        return filename


@click.command()
@click.argument("peptides")
@click.option(
    "--peptide-ksize",
    default=None,
    type=int,
    help="K-mer size of the peptide sequence to use. Defaults for"
    " different molecules are, "
    f"protein: {constants_index.DEFAULT_PROTEIN_KSIZE}"
    f", dayhoff: {constants_index.DEFAULT_DAYHOFF_KSIZE},"
    f" hydrophobic-polar: {constants_index.DEFAULT_HP_KSIZE}",
)
@click.option(
    "--alphabet",
    "--molecule",
    default="protein",
    help="The type of amino acid alphabet/encoding to use. Default "
    "is 'protein', but 'dayhoff' or 'hydrophobic-polar' can "
    "be used",
)
@click.option(
    "--save-as",
    default=None,
    help="If provided, save peptide bloom filter as this filename. "
    "Otherwise, add ksize and alphabet name to input filename.",
)
@click.option(
    "--tablesize",
    type=constants_index.BASED_INT,
    default="1e8",
    help="Size of the bloom filter table to use",
)
@click.option(
    "--n-tables",
    type=int,
    default=constants_index.DEFAULT_N_TABLES,
    help="Size of the bloom filter table to use",
)
@click.option(
    "--index-from-dir",
    is_flag=True,
    default=False,
    help="Build peptide bloom filter from a directory of peptide fasta files",
)
def cli(
    peptides,
    peptide_ksize=None,
    alphabet="protein",
    save_as=None,
    tablesize=constants_index.DEFAULT_MAX_TABLESIZE,
    n_tables=constants_index.DEFAULT_N_TABLES,
    index_from_dir=False,
):
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
    index_from_dir : bool
        peptides is a directory containing peptide sequence files
    long_reads
    verbose

    \b
    Returns
    -------

    """
    # \b above prevents rewrapping of paragraph

    index_dir = ""
    if index_from_dir:
        index_dir = peptides
        peptides = (os.path.join(index_dir, p) for p in os.listdir(index_dir))

    peptide_ksize = get_peptide_ksize(alphabet, peptide_ksize)
    peptide_bloom_filter = make_peptide_bloom_filter(
        peptides,
        peptide_ksize,
        alphabet,
        n_tables=n_tables,
        tablesize=tablesize,
        index_dir=index_dir,
    )
    logger.info("\tDone!")

    save_peptide_bloom_filter = save_as if save_as is not None else True
    maybe_save_peptide_bloom_filter(
        peptides,
        peptide_bloom_filter,
        alphabet,
        save_peptide_bloom_filter=save_peptide_bloom_filter,
        index_dir=index_dir,
    )


def get_peptide_ksize(molecule, peptide_ksize=None):
    if peptide_ksize is None:
        try:
            peptide_ksize = BEST_KSIZES[molecule]
        except KeyError:
            raise ValueError(
                f"{molecule} does not have a default k-mer size! "
                f"Only 'protein', 'hydrophobic-polar', or"
                f" 'dayhoff' have a default protein ksize"
            )
    return peptide_ksize
