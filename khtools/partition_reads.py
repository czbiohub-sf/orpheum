"""
partition_reads.py

Partition reads into coding, noncoding, and low-complexity bins
"""
import os
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

def write_fasta(file_handle, description, sequence):
    file_handle.write(f">{description}\n{sequence}\n")


def open_and_announce(filename, seqtype, quiet=False):
    if not quiet:
        print(f"Writing {seqtype} to {filename}")
    return open(filename, 'w')


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
        peptide_graph = make_peptide_bloom_filter(peptides, peptide_ksize,
                                                  molecule=molecule)
    else:
        peptide_graph = Nodegraph.load(peptides)
        assert peptide_ksize == peptide_graph.ksize
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


def three_frame_translation_no_stops(seq, debug=False, sign=1):
    """Remove translations with stop codons and keep track of the reading frame"""
    return {sign*i: t for i, t in enumerate(three_frame_translation(seq, debug))
            if '*' not in t}


def six_frame_translation_no_stops(seq, debug=False):
    forward_translations = three_frame_translation_no_stops(seq, debug)

    # Sign=-1 sets the reading frames as negative to make it obvious they are
    # from the reverse strand
    reverse_translations = three_frame_translation_no_stops(
        seq.reverse_complement(), debug, sign=-1)
    forward_translations.update(reverse_translations)
    return forward_translations


def score_single_translation(translation, peptide_graph, peptide_ksize,
                             molecule='protein', jaccard_threshold=0.9,
                             verbose=True, description=None,
                             translation_frame=None,
                             peptide_file_handle=None):
    translation = encode_peptide(translation, molecule)

    if len(translation) <= peptide_ksize:
        return 0, 0

    kmers = list(set(kmerize(str(translation), peptide_ksize)))
    hashes = [hash_murmur(kmer) for kmer in kmers]
    n_kmers = len(kmers)
    n_kmers_in_peptide_db = sum(1 for h in hashes if
                                peptide_graph.get(h) > 0)
    if verbose > 1:
        print(f"\ttranslation: \t{translation}")
        print("\tkmers:", ' '.join(kmers))

    kmers_in_peptide_db = {(k, h): peptide_graph.get(h) for k, h in
                           zip(kmers, hashes)}
    if verbose > 1:
        # Print keys (kmers) only
        print(f"\tK-mers in peptide database:")
        pprint(kmers_in_peptide_db)

    fraction_in_peptide_db = n_kmers_in_peptide_db / n_kmers

    if fraction_in_peptide_db > jaccard_threshold:
        if verbose:
            print(f"\t{translation} is above {jaccard_threshold}")
        if peptide_file_handle is not None:
            seqname = f'{description} translation_frame: {translation_frame}'
            write_fasta(peptide_file_handle, seqname, translation)

    return fraction_in_peptide_db, n_kmers


def compute_low_complexity(sequence, ksize):
    """Check if sequence is low complexity, i.e. low entropy, mostly repetitive"""
    kmers = kmerize(sequence, ksize)
    n_kmers = len(kmers)
    n_possible_kmers_on_sequence = len(sequence) - ksize + 1
    min_kmer_entropy = n_possible_kmers_on_sequence / 2
    if n_kmers < min_kmer_entropy:
        return True, n_kmers
    return False, n_kmers


def score_single_sequence(sequence, peptide_graph, peptide_ksize,
                          molecule='protein', verbose=True,
                          jaccard_threshold=0.9,
                          description=None,
                          peptide_file_handle=None,
                          noncoding_file_handle=None,
                          low_complexity_peptide_file_handle=None):
    # Convert to BioPython sequence object for translation
    seq = Seq(sequence)

    # Convert to BioPython sequence object for translation
    translations = six_frame_translation_no_stops(seq)
    # For all translations, use the one with the maximum number of k-mers
    # in the databse
    max_n_kmers = 0
    max_fraction_in_peptide_db = 0
    for frame, translation in translations.items():
        translation = encode_peptide(translation, molecule)

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            # Ignore Biopython warning of seq objects being strings now
            low_complexity, n_kmers = compute_low_complexity(translation,
                                                             peptide_ksize)
        if low_complexity:
            if peptide_file_handle is not None:
                seqname = f'{description} translation_frame: {translation_frame}'
                write_fasta(low_complexity_peptide_file_handle, seqname,
                            translation)
            return -1, n_kmers

        fraction_in_peptide_db, n_kmers = score_single_translation(
            translation, peptide_graph, peptide_ksize, molecule=molecule,
            verbose=verbose, description=description, translation_frame=frame,
            peptide_file_handle=peptide_file_handle)

        # Save the highest jaccard
        max_fraction_in_peptide_db = max(max_fraction_in_peptide_db,
                                         fraction_in_peptide_db)

        if max_fraction_in_peptide_db == fraction_in_peptide_db:
            # Update n_kmers if this is the best translation frame
            max_n_kmers = n_kmers

    if max_fraction_in_peptide_db <= jaccard_threshold:
        if noncoding_file_handle is not None:
            write_fasta(noncoding_file_handle, description, sequence)
    return max_fraction_in_peptide_db, max_n_kmers



def score_reads(reads, peptide_graph, peptide_ksize, jaccard_threshold=0.9,
                molecule='protein', verbose=False, prefix=None):
    scoring_lines = []
    nucleotide_ksize = 3*peptide_ksize

    if prefix is not None:
        noncoding_file_handle = open_and_announce(
            f"{prefix}.noncoding_nucleotides.fasta", "Noncoding nucleotides")
        peptide_file_handle = open_and_announce(
            f"{prefix}.coding_peptides.fasta",
             "all valid protein-coding translation frames")
        low_complexity_file_handle = open_and_announce(
            f"{prefix}.low_complexity_nucleotides.fasta",
             "low complexity (low entropy) nucleotides")
        low_complexity_peptide_file_handle = open_and_announce(
            f"{prefix}.low_complexity_peptides.fasta",
             "low complexity (low entropy) peptides")
    else:
        noncoding_file_handle, peptide_file_handle = None, None
        low_complexity_file_handle = None
        low_complexity_peptide_file_handle = None

    for record in tqdm(screed.open(reads)):
        description = record['name']
        sequence = record['sequence']
        if verbose:
            print(description)

        # Check if nucleotide sequence is low complexity
        is_low_complexity, n_kmers = compute_low_complexity(sequence,
                                                            nucleotide_ksize)
        if is_low_complexity:
            scoring_lines.append(
                [description, -1, n_kmers, 'low complexity'])
            if prefix is not None:
                write_fasta(low_complexity_file_handle, description, sequence)
            continue

        jaccard, n_kmers = score_single_sequence(
            sequence, peptide_graph, peptide_ksize, molecule, verbose,
            jaccard_threshold=jaccard_threshold,
            description=description,
            peptide_file_handle=peptide_file_handle,
            noncoding_file_handle=noncoding_file_handle,
            low_complexity_peptide_file_handle=low_complexity_peptide_file_handle)

        if jaccard > jaccard_threshold:
            line = [description, jaccard, n_kmers, 'coding']
        else:
            line = [description, jaccard, n_kmers, 'non-coding']
        if verbose > 1:
            # pprint(n_kmers)
            print(f"Jaccard: {jaccard}, n_kmers = {n_kmers}")
        scoring_lines.append(line)

    if prefix:
        noncoding_file_handle.close()
        peptide_file_handle.close()
        low_complexity_file_handle.close()
        low_complexity_peptide_file_handle.close()

    scoring_df = pd.DataFrame(scoring_lines,
                                 columns=['read_id',
                                          'jaccard_in_peptide_db',
                                          'n_kmers',
                                          'classification'])
    return scoring_df


def maybe_save_peptide_bloom_filter(peptides, peptide_graph,
                                    save_peptide_bloom_filter):
    if save_peptide_bloom_filter:
        if isinstance(save_peptide_bloom_filter, str):
            peptide_graph.save(save_peptide_bloom_filter)
        else:
            filename = os.path.splitext(peptides)[0] + '.nodegraph'
            print(f"Writing peptide bloom filter to {filename}")
            peptide_graph.save(filename)

@click.command()
@click.argument('reads')
@click.argument('peptides')
@click.option('--peptide-ksize', default=7,
                help="K-mer size of the peptide sequence to use. Default: 7")
@click.option("--save-peptide-bloom-filter", is_flag=True,
              help="If specified, save the peptide bloom filter. "
                   "Default filename is the name of the")
@click.option('--peptides-are-bloom-filter', is_flag=True,
              help="Peptide file is already a bloom filter")
@click.option('--jaccard-threshold', default=0.9,
              help="Minimum fraction of peptide k-mers from read in the "
                   "peptide database for this read to be called a "
                   "'coding read'")
@click.option('--molecule', default='protein',
              help="The type of amino acid encoding to use. Default is "
                   "'protein', but 'dayhoff' or 'hydrophobic-polar' can be "
                   "used")
@click.option('--csv', default=False,
               help='Name of csv file to write with all sequence reads and '
                    'their coding scores')
@click.option("--long-reads", is_flag=True,
              help="If set, then only considers reading frames starting with "
                   "start codon (ATG) and ending in a stop codon "
                   "(TAG, TAA, TGA)")
@click.option("--verbose", is_flag=True,
              help="Print more output")
@click.option("--debug", is_flag=True,
                  help="Print developer debugger output, including warnings")
def cli(reads, peptides, peptide_ksize=7, save_peptide_bloom_filter=True,
        peptides_are_bloom_filter=False, jaccard_threshold=0.9,
        molecule='protein', csv=False, long_reads=False, verbose=False,
        debug=False):
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
    # \b above prevents rewrapping of paragraph
    click.echo(f"Creating peptide bloom filter with file: {peptides} using " \
              f"ksize: {ksize} and molecule: {molecule} ...")
    peptide_graph = maybe_make_peptide_bloom_filter(peptides, peptide_ksize,
                                                    molecule,
                                                    peptides_are_bloom_filter)
    click.echo("\tDone!")

    prefix = os.path.splitext(reads)[0]
    coding_scores = score_reads(reads, peptide_graph, peptide_ksize,
                                jaccard_threshold, molecule, verbose,
                                prefix=prefix)
    if csv:
        coding_scores.to_csv(csv)

    maybe_save_peptide_bloom_filter(peptides, peptide_graph,
                                    save_peptide_bloom_filter)

