"""
partition_reads.py

Partition reads into coding, noncoding, and low-complexity bins
"""
from pprint import pprint
import warnings

from Bio.Seq import Seq
import click
from khmer import Nodegraph
import screed
from sourmash._minhash import hash_murmur
import pandas as pd
from khtools.sequence_encodings import encode_peptide
from khtools.compare_kmer_content import kmerize
from tqdm import tqdm


# Import modified 'os' module with LC_LANG set so click doesn't complain.
# The '# noqa: F401' line prevents the linter from complaining about the unused
# import.


DEFAULT_K = 32
DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = 1e10
DEFAULT_N_THREADS = 1
DEFAULT_SEED = 42



def make_peptide_bloom_filter(peptide_fasta, peptide_ksize, n_tables=4,
                              molecule='protein',
                              tablesize=DEFAULT_MAX_TABLESIZE):
    """Create a bloom filter out of peptide seuqences"""
    peptide_graph = Nodegraph(peptide_ksize, tablesize, n_tables=n_tables)

    for record in screed.open(peptide_fasta):
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


def three_frame_translation(seq, debug=False):
    if debug:
        warning_filter = 'default'
    else:
        warning_filter = 'ignore'

    with warnings.catch_warnings():
        warnings.simplefilter(warning_filter)
        for frame in range(3):
            translation = seq[frame:].translate()
            yield translation


def three_frame_translation_no_stops(seq, debug=False):
    return [t for t in three_frame_translation(seq, debug)
            if '*' not in t]


def six_frame_translation_no_stops(seq, debug=False):
    forward_translations = three_frame_translation_no_stops(seq, debug)
    reverse_translations = three_frame_translation_no_stops(
        seq.reverse_complement(), debug)
    return forward_translations + reverse_translations


def score_single_translation(translation, peptide_graph, peptide_ksize,
                             molecule='protein', jaccard_threshold=0.9,
                             verbose=True):
    translation = encode_peptide(translation, molecule)

    if len(translation) < peptide_ksize:
        return 0, 0

    kmers = list(set(kmerize(str(translation), peptide_ksize)))
    hashes = [hash_murmur(kmer) for kmer in kmers]
    n_kmers = len(kmers)
    n_kmers_in_peptide_db = sum(1 for h in hashes if
                                peptide_graph.get(h) > 0)
    if verbose:
        print(f"\ttranslation: \t{translation}")
        print("\tkmers:", ' '.join(kmers))

    if n_kmers < (len(translation) - peptide_ksize + 1) / 2:
        return -1
        if verbose:
            print(f'Low complexity sequence!!! n_kmers < (len(read.seq) - ksize + 1)/2  --> {n_kmers} < {(len(record.seq) - ksize + 1)/2}')
            print(record.description)
            print(record.seq)

    kmers_in_peptide_db = {(k, h): peptide_graph.get(h) for k, h in
                           zip(kmers, hashes)}
    if verbose:
        # Print keys (kmers) only
        print(f"\tK-mers in peptide database:")
        pprint(kmers_in_peptide_db)

    fraction_in_peptide_db = n_kmers_in_peptide_db / n_kmers

    if verbose and fraction_in_peptide_db > jaccard_threshold:
        print(f"{translation} is above {jaccard_threshold}")
    return fraction_in_peptide_db, n_kmers


def is_low_complexity(sequence, ksize):
    """CHeck if seqauence is low complexity, i.e. low entropy, mostly repetitive"""
    kmers = kmerize(sequence, ksize)
    n_kmers = len(kmers)
    n_possible_kmers_on_sequence = len(sequence) - ksize + 1
    min_kmer_entropy = n_possible_kmers_on_sequence / 2
    if n_kmers < min_kmer_entropy:
        return True, n_kmers
    return False, n_kmers


def score_single_sequence(sequence, peptide_graph, peptide_ksize,
                          molecule='protein', verbose=True):
    nucleotide_ksize = 3 * peptide_ksize
    # Check if nucleotide sequence is low complexity
    low_complexity, n_kmers = is_low_complexity(sequence,
                                                nucleotide_ksize)
    if low_complexity:
        return -1, n_kmers

    # Convert to BioPython sequence object for translation
    seq = Seq(sequence)

    # Convert to BioPython sequence object for translation
    translations = six_frame_translation_no_stops(seq)
    # For all translations, use the one with the maximum number of k-mers
    # in the databse
    max_n_kmers = 0
    max_fraction_in_peptide_db = 0
    for translation in translations:
        translation = encode_peptide(translation, molecule)
        fraction_in_peptide_db, n_kmers = score_single_translation(
            translation, peptide_graph, peptide_ksize, molecule=molecule,
            verbose=verbose)

        # Save the highest jaccard
        max_fraction_in_peptide_db = max(max_fraction_in_peptide_db,
                                         fraction_in_peptide_db)

        if max_fraction_in_peptide_db == fraction_in_peptide_db:
            # Update n_kmers if this is the best translation frame
            max_n_kmers = n_kmers
    return max_fraction_in_peptide_db, max_n_kmers



def score_reads(reads, peptide_graph, peptide_ksize, jaccard_threshold=0.9,
                molecule='protein', verbose=False):
    scoring_lines = []
    nucleotide_ksize = 3*peptide_ksize

    for record in tqdm(screed.open(reads)):
        description = record['name']
        sequence = record['sequence']

        # Check if nucleotide sequence is low complexity
        low_complexity, n_kmers = is_low_complexity(sequence,
                                                    nucleotide_ksize)
        if low_complexity:
            scoring_lines.append(
                [description, -1, n_kmers, 'low complexity'])
            continue

        jaccard, n_kmers = score_single_sequence(sequence, peptide_graph,
                                                 peptide_ksize, molecule, verbose)

        if jaccard > jaccard_threshold:
            line = [description, jaccard, n_kmers, 'coding']
        else:
            line = [description, jaccard, n_kmers, 'non-coding']
        if verbose:
            # pprint(n_kmers)
            print("Jaccard: {jaccard}, n_kmers = {n_kmers}")
        scoring_lines.append(line)

    scoring_df = pd.DataFrame(scoring_lines,
                                 columns=['read_id',
                                          'jaccard_in_peptide_db',
                                          'n_kmers',
                                          'classification'])
    return scoring_df


@click.command()
@click.argument('reads')
@click.argument('peptides')
@click.option('--peptide-ksize', default=7,
                help="K-mer size of the peptide sequence to use")
@click.option('--jaccard-threshold', default=0.9,
              help="Minimum fraction of peptide k-mers from read in the "
                   "peptide database for this read to be called a "
                   "'coding read'")
@click.option('--molecule', default='protein',
              help="The type of amino acid encoding to use. Default is "
                   "'protein', but 'dayhoff' or 'hydrophobic-polar' can be "
                   "used")
@click.option("--long-reads", is_flag=True,
              help="If set, then only considers reading frames starting with "
                   "start codon (ATG) and ending in a stop codon "
                   "(TAG, TAA, TGA)")
@click.option("--verbose", is_flag=True,
              help="Print more output")
@click.option("--debug", is_flag=True,
                  help="Print developer debugger output, including warnings")
def cli(reads, peptides, peptide_ksize, jaccard_threshold=0.9,
        molecule='protein', long_reads=False, verbose=False):
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
    coding_scores = score_reads(reads, peptide_graph, peptide_ksize,
                                jaccard_threshold, molecule, verbose)



