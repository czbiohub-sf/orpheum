"""
partition_reads.py

Partition reads into coding, noncoding, and low-complexity bins
"""
import os
from pprint import pprint
import warnings

from Bio.Seq import Seq
import click
import screed
from sourmash._minhash import hash_murmur
import pandas as pd
from khtools.sequence_encodings import encode_peptide
from khtools.compare_kmer_content import kmerize
from khtools.bloom_filter import (maybe_make_peptide_bloom_filter,
                                  maybe_save_peptide_bloom_filter)
from tqdm import tqdm


# Import modified 'os' module with LC_LANG set so click doesn't complain.
# The '# noqa: F401' line prevents the linter from complaining about the unused
# import.
DEFAULT_JACCARD_THRESHOLD = 0.5


def write_fasta(file_handle, description, sequence):
    file_handle.write(f">{description}\n{sequence}\n")


def open_and_announce(filename, seqtype, quiet=False):
    if not quiet:
        click.echo(f"Writing {seqtype} to {filename}", err=True)
    return open(filename, 'w')


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
    return {sign*(i+1): t for i, t in
            enumerate(three_frame_translation(seq, debug)) if '*' not in t}


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
        click.echo(f"\ttranslation: \t{translation}", err=True)
        click.echo("\tkmers:", ' '.join(kmers), err=True)

    if verbose > 1:
        kmers_in_peptide_db = {(k, h): peptide_graph.get(h) for k, h in
                               zip(kmers, hashes)}
        # Print keys (kmers) only
        click.echo(f"\tK-mers in peptide database:", err=True)
        click.echo(kmers_in_peptide_db, err=True)

    fraction_in_peptide_db = n_kmers_in_peptide_db / n_kmers

    if fraction_in_peptide_db > jaccard_threshold:
        if verbose:
            click.echo(f"\t{translation} is above {jaccard_threshold}",
                        err=True)
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
                          translation_frame=None,
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
        # Convert back to string
        translation = str(translation)

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            # Ignore Biopython warning of seq objects being strings now
            low_complexity, n_kmers = compute_low_complexity(translation,
                                                             peptide_ksize)

        # Maybe reencode to dayhoff/hp space
        encoded = encode_peptide(translation, molecule)
        if low_complexity:
            if low_complexity_peptide_file_handle is not None:
                seqname = f'{description} translation_frame: {frame}'
                write_fasta(low_complexity_peptide_file_handle, seqname,
                            translation)
            return -2, n_kmers

        fraction_in_peptide_db, n_kmers = score_single_translation(
            encoded, peptide_graph, peptide_ksize, molecule=molecule,
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



def score_reads(reads, peptide_graph, peptide_ksize,
                jaccard_threshold=DEFAULT_JACCARD_THRESHOLD,
                molecule='protein', verbose=False, prefix=None):
    scoring_lines = []
    nucleotide_ksize = 3*peptide_ksize

    if prefix is not None:
        noncoding_file_handle = open_and_announce(
            f"{prefix}.noncoding_nucleotides.fasta", "noncoding nucleotides")
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
                [description, -1, n_kmers, 'low complexity nucleotide'])
            if low_complexity_file_handle is not None:
                write_fasta(low_complexity_file_handle, description, sequence)
            continue

        jaccard, n_kmers = score_single_sequence(
            sequence, peptide_graph, peptide_ksize, molecule, verbose,
            jaccard_threshold=jaccard_threshold,
            description=description,
            peptide_file_handle=peptide_file_handle,
            noncoding_file_handle=noncoding_file_handle,
            low_complexity_peptide_file_handle=low_complexity_peptide_file_handle)

        if jaccard == -2:
            line = [description, jaccard, n_kmers, 'low complexity peptide']
        elif jaccard > jaccard_threshold:
            line = [description, jaccard, n_kmers, 'coding']
        else:
            line = [description, jaccard, n_kmers, 'non-coding']
        if verbose > 1:
            # pprint(n_kmers)
            click.echo(f"Jaccard: {jaccard}, n_kmers = {n_kmers}", err=True)
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


@click.command()
@click.argument('peptides', nargs=1)
@click.argument('reads', nargs=-1)
@click.option('--peptide-ksize', default=7,
                help="K-mer size of the peptide sequence to use. Default: 7")
@click.option("--save-peptide-bloom-filter", is_flag=True, default=False,
              help="If specified, save the peptide bloom filter. "
                   "Default filename is the name of the")
@click.option('--peptides-are-bloom-filter', is_flag=True, default=False,
              help="Peptide file is already a bloom filter")
@click.option('--jaccard-threshold', default=DEFAULT_JACCARD_THRESHOLD,
              help="Minimum fraction of peptide k-mers from read in the "
                   "peptide database for this read to be called a " +
                   f"'coding read'. Default: {DEFAULT_JACCARD_THRESHOLD")
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
def cli(peptides, reads, peptide_ksize=7, save_peptide_bloom_filter=True,
        peptides_are_bloom_filter=False,
        jaccard_threshold=DEFAULT_JACCARD_THRESHOLD,
        molecule='protein', csv=False, long_reads=False, verbose=False,
        debug=False):
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
    peptide_graph = maybe_make_peptide_bloom_filter(peptides, peptide_ksize,
                                                    molecule,
                                                    peptides_are_bloom_filter)
    click.echo("\tDone!")

    peptides_basename = os.path.basename(peptides)
    suffix = f"__{peptides_basename}__molecule-{molecule}__" \
             f"ksize-{peptide_ksize}"

    if not peptides_are_bloom_filter:
        maybe_save_peptide_bloom_filter(peptides, peptide_graph,
                                        molecule, peptide_ksize,
                                        save_peptide_bloom_filter)

    dfs = []
    for reads_file in reads:
        prefix = os.path.splitext(reads_file)[0] + suffix
        df = score_reads(reads_file, peptide_graph, peptide_ksize,
                                    jaccard_threshold, molecule, verbose,
                                    prefix=prefix)
        df[filename] = reads_file
        dfs.append(df)

    coding_scores = pd.concat(dfs, ignore_index=True)

    if csv:
        click.echo(f"Writing coding scores of reads to {csv}", err=True)
        coding_scores.to_csv(csv, index=False)


if __name__ == '__main__':
    cli()
