import math
import os

import click
import khmer
import screed
from sourmash._minhash import hash_murmur
from tqdm import tqdm

from sencha.compare_kmer_content import kmerize
from sencha.sequence_encodings import encode_peptide, ALPHABET_SIZES
from sencha.constants_index import BEST_KSIZES
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


def make_peptide_index(
    peptide_fasta,
    peptide_ksize,
    molecule,
    n_tables=constants_index.DEFAULT_N_TABLES,
    tablesize=constants_index.DEFAULT_MAX_TABLESIZE,
    max_observed_fraction=constants_index.MAX_FRACTION_OBSERVED_TO_THEORETICAL_KMERS,
):
    """Create a bloom filter out of peptide sequences"""
    peptide_index = khmer.Nodegraph(peptide_ksize, tablesize, n_tables=n_tables)

    with screed.open(peptide_fasta) as records:
        for record in tqdm(records):
            sequence = encode_peptide(record["sequence"], molecule)
            if len(sequence) >= peptide_ksize:
                # Skip sequences with any stop codons
                if '*' in sequence:
                    continue
                kmers = kmerize(sequence, peptide_ksize)
                for kmer in kmers:
                    # Ignore the k-mer if there are any illegal characters
                    if any(x in kmer for x in constants_index.RESIDUES_TO_IGNORE):
                        logger.info(
                            f'{record["name"]} a k-mer ({kmer}) with an '
                            f"illegal character, skipping this k-mer"
                        )
                        continue
                    # Convert the k-mer into an integer
                    hashed = hash_murmur(kmer)

                    # .add can take the hashed integer so we can hash the
                    #  peptide kmer and add it directly
                    peptide_index.add(hashed)
            else:
                logger.info(
                    f'{record["name"]} sequence is shorter than the k-mer '
                    f"size {peptide_ksize}, skipping"
                )
    collisions = khmer.calc_expected_collisions(peptide_index, force=True)
    if collisions > constants_index.MAX_BF_FALSE_POSITIVES:
        raise ValueError(f"The false positive rate in the bloom filter index is "
                         f"{collisions}, which is greater than than the recommended "
                         f"maximum of {constants_index.MAX_BF_FALSE_POSITIVES:.1f}. "
                         f"The current table size is {tablesize:.1e}, please increase "
                         f"by an order of magnitude and rerun.")

    n_theoretical_kmers = ALPHABET_SIZES[molecule] ** peptide_ksize
    n_observed_kmers = peptide_index.n_unique_kmers()
    fraction_observed = n_observed_kmers / n_theoretical_kmers
    if fraction_observed > max_observed_fraction:
        raise ValueError(
            f"The number of observed length {peptide_ksize} k-mers compared to the "
            f"possible theoretical k-mers "
            f"is {n_observed_kmers} / {n_theoretical_kmers} = {fraction_observed:.2e} "
            f"which is greater than the maximum observed fraction threshold, "
            f"{max_observed_fraction:.2e}. "
            f"This doesn't leave enough 'negative space' of non-observed protein "
            f"k-mers for room for predicting true protein-coding sequence, which "
            f"relies on seeing which protein k-mers are *not* present in the observed "
            f"data. Please increase the k-mer size."
        )
    return peptide_index


def make_peptide_set(peptide_fasta, peptide_ksize, molecule):
    """Create a python set out of peptide sequence k-mers

    For comparing to the bloom filter in storage and performance
    """
    peptide_set = set([])

    with screed.open(peptide_fasta) as records:
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


def maybe_make_peptide_index(
    peptides,
    peptide_ksize,
    molecule,
    peptides_are_index,
    n_tables=constants_index.DEFAULT_N_TABLES,
    tablesize=constants_index.DEFAULT_MAX_TABLESIZE,
):
    if peptides_are_index:
        logger.info(
            f"Loading existing bloom filter from {peptides} and "
            f"making sure the ksizes match"
        )
        peptide_index = load_nodegraph(peptides)
        if peptide_ksize is not None:
            try:
                assert peptide_ksize == peptide_index.ksize()
            except AssertionError:
                raise ValueError(
                    f"Given peptide ksize ({peptide_ksize}) and "
                    f"ksize found in bloom filter "
                    f"({peptide_index.ksize()}) are not"
                    f"equal"
                )
    else:
        peptide_ksize = get_peptide_ksize(molecule, peptide_ksize)
        logger.info(
            f"Creating peptide bloom filter with file: {peptides}\n"
            f"Using ksize: {peptide_ksize} and alphabet: {molecule} "
            f"..."
        )
        peptide_index = make_peptide_index(
            peptides,
            peptide_ksize,
            molecule=molecule,
            n_tables=n_tables,
            tablesize=tablesize,
        )
    return peptide_index


def maybe_save_peptide_index(
    peptides, peptide_bloom_filter, molecule, save_peptide_index
):
    if save_peptide_index:
        ksize = peptide_bloom_filter.ksize()

        if isinstance(save_peptide_index, str):
            filename = save_peptide_index
            peptide_bloom_filter.save(save_peptide_index)
        else:
            suffix = f".alphabet-{molecule}_ksize-{ksize}.bloomfilter." f"nodegraph"
            filename = os.path.splitext(peptides)[0] + suffix

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
    "--max-observed-fraction",
    type=float,
    default=constants_index.MAX_FRACTION_OBSERVED_TO_THEORETICAL_KMERS,
    help="Maximum fraction of observed to theoretical k-mers to allow. Theoretical "
         "k-mers are computed as the (alphabet size)^(peptide_ksize). This number "
         "should be fairly small (e.g. 1e-4) to prevent false positive translation "
         "results.",
)
def cli(
    peptides,
    peptide_ksize=None,
    alphabet="protein",
    save_as=None,
    tablesize=constants_index.DEFAULT_MAX_TABLESIZE,
    n_tables=constants_index.DEFAULT_N_TABLES,
    max_observed_fraction=constants_index.MAX_FRACTION_OBSERVED_TO_THEORETICAL_KMERS,
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
    long_reads
    verbose

    \b
    Returns
    -------

    """
    # \b above prevents rewrapping of paragraph
    peptide_ksize = get_peptide_ksize(alphabet, peptide_ksize)
    peptide_index = make_peptide_index(
        peptides, peptide_ksize, alphabet, n_tables=n_tables, tablesize=tablesize,
        max_observed_fraction=max_observed_fraction
    )
    logger.info("\tDone!")

    save_peptide_index = save_as if save_as is not None else True
    maybe_save_peptide_index(
        peptides, peptide_index, alphabet, save_peptide_index=save_peptide_index,
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
