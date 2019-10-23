"""
extract_coding.py

Partition reads into coding, noncoding, and low-complexity bins
"""
import os
import sys
import warnings

from Bio.Seq import Seq
import click
import numpy as np
import pandas as pd
import screed
from sourmash._minhash import hash_murmur
from khtools.sequence_encodings import encode_peptide
from khtools.compare_kmer_content import kmerize
from khtools.bloom_filter import (maybe_make_peptide_bloom_filter,
                                  maybe_save_peptide_bloom_filter,
                                  get_peptide_ksize, DEFAULT_PROTEIN_KSIZE,
                                  DEFAULT_DAYHOFF_KSIZE, DEFAULT_HP_KSIZE)
from tqdm import tqdm


# Import modified 'os' module with LC_LANG set so click doesn't complain.
# The '# noqa: F401' line prevents the linter from complaining about the unused
# import.
DEFAULT_JACCARD_THRESHOLD = 0.5
SEQTYPE_TO_ANNOUNCEMENT = {"noncoding_nucleotide":
                               "nucleotide sequence from reads WITHOUT matches to "
                               "protein-coding peptides",
                           "coding_nucleotide":
                               "nucleotide sequence from reads WITH protein-coding translation"
                               " frame nucleotides",
                           "low_complexity_nucleotide":
                               "nucleotide sequence from low complexity (low entropy) reads",
                           "low_complexity_peptide": "peptide sequence from low "
                                                     "complexity (low entropy) translated"
                                                     " reads"
                           }
SCORING_DF_COLUMNS = ['read_id',
                      'jaccard_in_peptide_db',
                      'n_kmers',
                      'classification']


def write_fasta(file_handle, description, sequence):
    file_handle.write(f">{description}\n{sequence}\n")


def open_and_announce(filename, seqtype, quiet=False):
    if not quiet:
        announcement = SEQTYPE_TO_ANNOUNCEMENT[seqtype]
        click.echo(f"Writing {announcement} to {filename}", err=True)
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
    """Remove translations with stop codons & keep track of reading frame"""
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


def score_single_translation(translation, peptide_bloom_filter, peptide_ksize,
                             molecule='protein',
                             verbose=True):
    encoded = encode_peptide(translation, molecule)
    kmers = list(set(kmerize(str(encoded), peptide_ksize)))
    hashes = [hash_murmur(kmer) for kmer in kmers]
    n_kmers = len(kmers)
    n_kmers_in_peptide_db = sum(1 for h in hashes if peptide_bloom_filter.get(h) > 0)
    if verbose > 1:
        click.echo(f"\ttranslation: \t{encoded}", err=True)
        click.echo("\tkmers:", ' '.join(kmers), err=True)

    if verbose > 1:
        kmers_in_peptide_db = {(k, h): peptide_bloom_filter.get(h) for k, h in
                               zip(kmers, hashes)}
        # Print keys (kmers) only
        click.echo(f"\tK-mers in peptide database:", err=True)
        click.echo(kmers_in_peptide_db, err=True)

    fraction_in_peptide_db = n_kmers_in_peptide_db / n_kmers

    return fraction_in_peptide_db, n_kmers


def compute_low_complexity(sequence, ksize):
    """Check if sequence is low complexity, i.e. mostly repetitive"""
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        # Ignore Biopython warning of seq objects being strings now
        try:
            kmers = kmerize(sequence, ksize)
        except ValueError:
            # k-mer size is larger than sequence
            return True, 0
    n_kmers = len(kmers)
    n_possible_kmers_on_sequence = len(sequence) - ksize + 1
    min_kmer_entropy = n_possible_kmers_on_sequence / 2
    if n_kmers <= min_kmer_entropy:
        return True, n_kmers
    return False, n_kmers


def score_single_read(sequence, peptide_bloom_filter, peptide_ksize,
                      molecule='protein', verbose=True,
                      jaccard_threshold=0.9,
                      description=None,
                      noncoding_file_handle=None,
                      coding_nucleotide_file_handle=None,
                      low_complexity_peptide_file_handle=None):
    # Convert to BioPython sequence object for translation
    seq = Seq(sequence)

    # Convert to BioPython sequence object for translation
    translations = six_frame_translation_no_stops(seq)
    # For all translations, use the one with the maximum number of k-mers
    # in the databse
    max_n_kmers = 0
    max_fraction_in_peptide_db = 0
    if len(translations) == 0:
        return np.nan, np.nan, "No translation frames without stop codons"

    translations = {frame: translation for frame, translation
                    in translations.items()
                    if len(translation) > peptide_ksize}
    if len(translations) == 0:
        return np.nan, np.nan, "All translations shorter than peptide k-mer " \
                               "size + 1"

    for frame, translation in translations.items():
        # Convert back to string
        translation = str(translation)

        # Maybe reencode to dayhoff/hp space
        encoded = encode_peptide(translation, molecule)

        low_complexity, n_kmers = compute_low_complexity(encoded,
                                                         peptide_ksize)
        if low_complexity:
            if n_kmers > 0:
                maybe_write_fasta(description + f" translation_frame: {frame}",
                                  low_complexity_peptide_file_handle, translation)
                return np.nan, n_kmers, f"Low complexity peptide in {molecule}" \
                                        " encoding"
            else:
                return np.nan, np.nan, "Translated read length was smaller " \
                                       "than peptide k-mer size"

        fraction_in_peptide_db, n_kmers = score_single_translation(
            encoded, peptide_bloom_filter, peptide_ksize, molecule=molecule,
            verbose=verbose)

        # Save the highest jaccard
        max_fraction_in_peptide_db = max(max_fraction_in_peptide_db,
                                         fraction_in_peptide_db)

        if max_fraction_in_peptide_db == fraction_in_peptide_db:
            # Update n_kmers if this is the best translation frame
            max_n_kmers = n_kmers
        if fraction_in_peptide_db > jaccard_threshold:
            if verbose:
                click.echo(f"\t{translation} is above {jaccard_threshold}",
                           err=True)
            seqname = f'{description} translation_frame: {frame}'
            write_fasta(sys.stdout, seqname, translation)

    if max_fraction_in_peptide_db <= jaccard_threshold:
        maybe_write_fasta(description, noncoding_file_handle, sequence)
    return max_fraction_in_peptide_db, max_n_kmers, None


def maybe_write_fasta(description, file_handle, sequence):
    """Write fasta to file handle if it is not None"""
    if file_handle is not None:
        write_fasta(file_handle, description, sequence)


def score_reads(reads, peptide_bloom_filter, peptide_ksize,
                jaccard_threshold=DEFAULT_JACCARD_THRESHOLD,
                molecule='protein', verbose=False,
                coding_nucleotide_fasta=None,
                noncoding_nucleotide_fasta=None,
                low_complexity_nucleotide_fasta=None,
                low_complexity_peptide_fasta=None):
    """Assign a coding score to each read. Where the magic happens."""
    scoring_lines = []
    nucleotide_ksize = 3 * peptide_ksize

    fastas, file_handles = maybe_open_fastas(coding_nucleotide_fasta,
                                             low_complexity_nucleotide_fasta,
                                             low_complexity_peptide_fasta,
                                             noncoding_nucleotide_fasta)
    for record in tqdm(screed.open(reads)):
        description = record['name']
        sequence = record['sequence']
        if verbose:
            print(description)

        jaccard, n_kmers, special_case = maybe_score_single_read(
            description, fastas, file_handles, jaccard_threshold, molecule,
            nucleotide_ksize, peptide_bloom_filter, peptide_ksize, sequence,
            verbose)

        line = get_coding_score_line(description, jaccard, jaccard_threshold,
                                     n_kmers, special_case)
        scoring_lines.append(line)

    maybe_close_files(file_handles)

    # Concatenate all the lines into a single dataframe
    scoring_df = pd.DataFrame(scoring_lines, columns=SCORING_DF_COLUMNS)
    return scoring_df


def maybe_score_single_read(description, fastas, file_handles,
                            jaccard_threshold, molecule, nucleotide_ksize,
                            peptide_bloom_filter, peptide_ksize, sequence,
                            verbose):
    """Check if read is low complexity/too short, otherwise score it"""
    # Check if nucleotide sequence is low complexity
    is_low_complexity, n_kmers = compute_low_complexity(sequence,
                                                        nucleotide_ksize)
    if is_low_complexity:
        jaccard, n_kmers, special_case = too_short_or_low_complexity(
            description, fastas, n_kmers, sequence)
    else:
        jaccard, n_kmers, special_case = score_single_read(
            sequence, peptide_bloom_filter, peptide_ksize, molecule, verbose,
            jaccard_threshold=jaccard_threshold,
            description=description,
            noncoding_file_handle=file_handles['noncoding_nucleotide'],
            coding_nucleotide_file_handle=file_handles['coding_nucleotide'],
            low_complexity_peptide_file_handle=file_handles[
                'low_complexity_peptide'])

        if verbose > 1:
            click.echo(f"Jaccard: {jaccard}, n_kmers = {n_kmers}", err=True)
    return jaccard, n_kmers, special_case


def too_short_or_low_complexity(description, fastas, n_kmers,
                                sequence):
    if n_kmers > 0:
        jaccard = np.nan
        special_case = "Low complexity nucleotide"
        maybe_write_fasta(description, fastas['low_complexity_nucleotide'],
                          sequence)
    else:
        jaccard = np.nan
        n_kmers = 0
        special_case = 'Read length was shorter than 3 * preptide ' \
                       'k-mer size'
    return jaccard, n_kmers, special_case


def maybe_close_files(file_handles):
    for file_handle in file_handles.values():
        if file_handle is not None:
            file_handle.close()


def get_coding_score_line(description, jaccard, jaccard_threshold, n_kmers,
                          special_case):
    if special_case is not None:
        line = [description, jaccard, n_kmers, special_case]
    elif jaccard > jaccard_threshold:
        line = [description, jaccard, n_kmers, 'Coding']
    else:
        line = [description, jaccard, n_kmers, 'Non-coding']
    return line


def maybe_open_fastas(coding_nucleotide_fasta, low_complexity_nucleotide_fasta,
                      low_complexity_peptide_fasta,
                      noncoding_nucleotide_fasta):
    fastas = {"noncoding_nucleotide":
                  noncoding_nucleotide_fasta,
              "coding_nucleotide":
                  coding_nucleotide_fasta,
              "low_complexity_nucleotide":
                  low_complexity_nucleotide_fasta,
              "low_complexity_peptide":
                  low_complexity_peptide_fasta}
    file_handles = {}
    for seqtype, fasta in fastas.items():
        if fasta is not None:
            file_handles[seqtype] = open_and_announce(
                fasta, seqtype)
        else:
            file_handles[seqtype] = None
    return fastas, file_handles


@click.command()
@click.argument('peptides', nargs=1)
@click.argument('reads', nargs=-1)
@click.option('--peptide-ksize', default=None,
              help="K-mer size of the peptide sequence to use. Defaults for"
                     " different molecules are, "
                     f"protein: {DEFAULT_PROTEIN_KSIZE}"
                     f", dayhoff: {DEFAULT_DAYHOFF_KSIZE},"
                     f" hydrophobic-polar: {DEFAULT_HP_KSIZE}")
@click.option("--save-peptide-bloom-filter", is_flag=True, default=False,
              help="If specified, save the peptide bloom filter. "
                   "Default filename is the name of the fasta file plus a "
                   "suffix denoting the protein encoding and peptide ksize")
@click.option('--peptides-are-bloom-filter', is_flag=True, default=False,
              help="Peptide file is already a bloom filter")
@click.option('--jaccard-threshold', default=DEFAULT_JACCARD_THRESHOLD,
              help="Minimum fraction of peptide k-mers from read in the "
                   "peptide database for this read to be called a " +
                   f"'coding read'. Default: {DEFAULT_JACCARD_THRESHOLD}")
@click.option('--molecule', default='protein',
              help="The type of amino acid encoding to use. Default is "
                   "'protein', but 'dayhoff' or 'hydrophobic-polar' can be "
                   "used")
@click.option('--csv', default=False,
               help='Name of csv file to write with all sequence reads and '
                    'their coding scores')
@click.option("--coding-nucleotide-fasta",
              help="If specified, save the coding nucleotides to this file")
@click.option("--noncoding-nucleotide-fasta",
              help="If specified, save the noncoding nucleotides to this file")
@click.option("--low-complexity-nucleotide-fasta",
              help="If specified, save the low-complexity nucleotides to this file")
@click.option("--low-complexity-peptide-fasta",
              help="If specified, save the low-complexity peptides to this file")
@click.option("--long-reads", is_flag=True,
              help="If set, then only considers reading frames starting with "
                   "start codon (ATG) and ending in a stop codon "
                   "(TAG, TAA, TGA)")
@click.option("--verbose", is_flag=True,
              help="Print more output")
@click.option("--debug", is_flag=True,
                  help="Print developer debugger output, including warnings")
def cli(peptides, reads, peptide_ksize=None,
        save_peptide_bloom_filter=True,
        peptides_are_bloom_filter=False,
        jaccard_threshold=DEFAULT_JACCARD_THRESHOLD,
        molecule='protein', csv=False, long_reads=False,
        coding_nucleotide_fasta=None,
        noncoding_nucleotide_fasta=None,
        low_complexity_nucleotide_fasta=None,
        low_complexity_peptide_fasta=None,
        verbose=False,
        debug=False):
    """Writes coding peptides from reads to standard output

    Sane defaults for peptide_ksize for different peptide encodings:
    - with "protein" or "peptide" --> --peptide-ksize = 5-10
      7 is pretty universal but can go down to 5 for less species specificity
      and up to 10 to be very specific
    - with "dayhoff" --> --peptide-ksize = 10-15
    - with "hydrophobic-polar" or "hp" --> --peptide-ksize = 15-21
      15 is pretty good but can do up to 21

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

    if long_reads:
        raise NotImplementedError("Not implemented! ... yet :)")

    peptide_ksize = get_peptide_ksize(molecule, peptide_ksize)

    peptide_bloom_filter = maybe_make_peptide_bloom_filter(
        peptides, peptide_ksize, molecule, peptides_are_bloom_filter)
    click.echo("\tDone!")

    if not peptides_are_bloom_filter:
        maybe_save_peptide_bloom_filter(peptides, peptide_bloom_filter,
                                        molecule, peptide_ksize,
                                        save_peptide_bloom_filter)

    dfs = []
    for reads_file in reads:
        df = score_reads(
            reads_file, peptide_bloom_filter, peptide_ksize,
            jaccard_threshold, molecule, verbose,
            coding_nucleotide_fasta=coding_nucleotide_fasta,
            noncoding_nucleotide_fasta=noncoding_nucleotide_fasta,
            low_complexity_nucleotide_fasta=low_complexity_nucleotide_fasta,
            low_complexity_peptide_fasta=low_complexity_peptide_fasta)
        df[reads_file] = reads_file
        dfs.append(df)

    coding_scores = pd.concat(dfs, ignore_index=True)

    if csv:
        click.echo(f"Writing coding scores of reads to {csv}", err=True)
        coding_scores.to_csv(csv, index=False)


if __name__ == '__main__':
    cli()
