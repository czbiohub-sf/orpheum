"""
extract_coding.py

Partition reads into coding, noncoding, and low-complexity bins
"""
from collections import namedtuple
import itertools
import json
import sys
import warnings

from Bio.Seq import Seq
import click
import numpy as np
import pandas as pd
import screed
from sourmash._minhash import hash_murmur
from khtools.sequence_encodings import encode_peptide, ALPHABET_ALIASES
from khtools.compare_kmer_content import kmerize
from khtools.bloom_filter import (maybe_make_peptide_bloom_filter,
                                  maybe_save_peptide_bloom_filter,
                                  DEFAULT_PROTEIN_KSIZE,
                                  DEFAULT_DAYHOFF_KSIZE, DEFAULT_HP_KSIZE,
                                  DEFAULT_N_TABLES, DEFAULT_MAX_TABLESIZE,
                                  BASED_INT)
from tqdm import tqdm

# Import modified 'os' module with LC_LANG set so click doesn't complain.
# The '# noqa: F401' line prevents the linter from complaining about the unused
# import.
DEFAULT_JACCARD_THRESHOLD = 0.5
DEFAULT_HP_JACCARD_THRESHOLD = 0.8
SEQTYPE_TO_ANNOUNCEMENT = {
    "noncoding_nucleotide":
        "nucleotide sequence from reads WITHOUT matches to "
        "protein-coding peptides",
    "coding_nucleotide":
        "nucleotide sequence from reads WITH protein-coding translation"
        " frame nucleotides",
    "low_complexity_nucleotide":
        "nucleotide sequence from low complexity (low entropy) reads",
    "low_complexity_peptide":
        "peptide sequence from low "
        "complexity (low entropy) translated"
        " reads"
}
SCORING_DF_COLUMNS = [
    'read_id', 'jaccard_in_peptide_db', 'n_kmers', 'classification'
]

LOW_COMPLEXITY_PER_ALIAS = [
    list(
        (alias,
         f"Low complexity peptide in {alphabet} alphabet")
        for alias in aliases) for alphabet,
    aliases in ALPHABET_ALIASES.items()]
LOW_COMPLEXITY_CATEGORIES = dict(
    list(itertools.chain(*LOW_COMPLEXITY_PER_ALIAS)))


PROTEIN_CODING_CATEGORIES = {
    "too_short_peptide":
        "All translations shorter than peptide k-mer size + 1",
    "stop_codons": "All translation frames have stop codons",
    "coding": "Coding",
    'non_coding': 'Non-coding',
    'low_complexity_nucleotide': "Low complexity nucleotide",
    'too_short_nucleotide':
        'Read length was shorter than 3 * peptide k-mer size'
}

# Some mild object-orientation
SingleTranslationScore = namedtuple("SingleTranslationScore",
                                    ['fraction_kmers_in_peptide_db',
                                     'n_total_kmers'])
LowComplexityScore = namedtuple("LowComplexityScore",
                                ['is_low_complexity', 'n_total_kmers'])

SingleReadScore = namedtuple(
    "SingleReadScore",
    ['max_fraction_kmers_in_peptide_db_across_six_frame_translations',
     'max_n_kmers', 'special_case'])
PercentagesCounts = namedtuple("PercentagesCounts", ['percentages', 'counts'])
OutputFastasHandles = namedtuple('OutputFastasHandles',
                                 ['fastas', 'file_handles'])


def validate_jaccard(ctx, param, value):
    """Ensure Jaccard threshold is between 0 and 1"""
    if value is None:
        return value
    try:
        jaccard = float(value)
        assert jaccard <= 1
        assert jaccard >= 0
        return jaccard
    except (ValueError, AssertionError):
        raise click.BadParameter(f'--jaccard-threshold needs to be a number'
                                 f' between 0 and 1, but {value} was provided')


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
    return {
        sign * (i + 1): t
        for i, t in enumerate(three_frame_translation(seq, debug))
        if '*' not in t
    }


def six_frame_translation_no_stops(seq, debug=False):
    forward_translations = three_frame_translation_no_stops(seq, debug)

    # Sign=-1 sets the reading frames as negative to make it obvious they are
    # from the reverse strand
    reverse_translations = three_frame_translation_no_stops(
        seq.reverse_complement(), debug, sign=-1)
    forward_translations.update(reverse_translations)
    return forward_translations


def score_single_translation(translation,
                             peptide_bloom_filter,
                             peptide_ksize,
                             molecule='protein',
                             verbose=True):
    encoded = encode_peptide(translation, molecule)
    kmers = list(set(kmerize(str(encoded), peptide_ksize)))
    hashes = [hash_murmur(kmer) for kmer in kmers]
    n_kmers = len(kmers)
    n_kmers_in_peptide_db = sum(1 for h in hashes
                                if peptide_bloom_filter.get(h) > 0)
    if verbose > 1:
        click.echo(f"\ttranslation: \t{encoded}", err=True)
        click.echo("\tkmers:", ' '.join(kmers), err=True)

    if verbose > 1:
        kmers_in_peptide_db = {(k, h): peptide_bloom_filter.get(h)
                               for k, h in zip(kmers, hashes)}
        # Print keys (kmers) only
        click.echo(f"\tK-mers in peptide database:", err=True)
        click.echo(kmers_in_peptide_db, err=True)

    fraction_in_peptide_db = n_kmers_in_peptide_db / n_kmers
    return SingleTranslationScore(fraction_in_peptide_db, n_kmers)


def evaluate_is_fastp_low_complexity(seq, complexity_threshold=0.3):
    """Use fastp's definition of complexity

    By this definition, low complexity sequence is defined by consecutive runs
    of same base in a row, e.g.
    CCCCCCCCCACCACCACCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCACACACCCCCAACAC
    is low complexity. The threshold is 0.3 as used in the fastp prpject:
    https://github.com/OpenGene/fastp

    Parameters
    ----------
    seq : str
        Sequence to compute complexity on
    complexity_threshold : float, defaault 0.3
        Value between 0 and 1. The default is 0.3 because that is the default
        in the command line program fastp

    Returns
    -------
    is_low_complexity : bool
        Whether or not the sequence passes the complexity threshold
    """
    complexity = compute_fastp_complexity(seq)
    return complexity < complexity_threshold


def compute_fastp_complexity(seq):
    n_different_consecutively = sum(1 for i in range(len(seq) - 1)
                                    if seq[i] != seq[i + 1])
    complexity = n_different_consecutively / len(seq)
    return complexity


def evaluate_is_kmer_low_complexity(sequence, ksize):
    """Check if sequence is low complexity, i.e. mostly repetitive

    By this definition, the sequence is not complex if its number of unique
    k-mers is smaller than half the number of expected k-mers
    """
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        # Ignore Biopython warning of seq objects being strings now
        try:
            kmers = kmerize(sequence, ksize)
        except ValueError:
            # k-mer size is larger than sequence
            return LowComplexityScore(None, None)
    n_kmers = len(kmers)
    n_possible_kmers_on_sequence = len(sequence) - ksize + 1
    min_kmer_entropy = n_possible_kmers_on_sequence / 2
    is_low_complexity = n_kmers <= min_kmer_entropy
    return LowComplexityScore(is_low_complexity, n_kmers)


def score_single_read(sequence,
                      peptide_bloom_filter,
                      peptide_ksize,
                      alphabet='protein',
                      verbose=True,
                      jaccard_threshold=0.9,
                      description=None,
                      noncoding_file_handle=None,
                      coding_nucleotide_file_handle=None,
                      low_complexity_peptide_file_handle=None):
    """Predict whether a nucleotide sequence could be protein-coding

    Parameters
    ----------
    sequence : str
        Nucleotide sequence to predict on
    peptide_bloom_filter : khmer.Nodegraph
        Database of known peptide k-mers from a well-studied organism, e.g.
        human protein-coding sequences. Must have been built on peptides using
        the same k-mer size and molecular encoding as specified here, otherwise
        the results will make no sense
    peptide_ksize : int
        Length of the peptide words in sequence. Must match the k-mer size used
        for the peptide_bloom_filter otherwise nothing will match, or only
        false positives will match.
    alphabet : str
        One of "protein"|"peptide", "dayhoff", or "hydrophobic-polar"|"hp" to
        encode the protein-coding space. Where "protein"|"peptide" is the
        original 20-letter amino acid encoding, Dayhoff ("dayhoff") is a lossy
        6-letter encoding that categorizes the amino acids into:
            1. Cysteine,
            2. Small (A, G, P, S, T)
            3. Acid and Amide (D, E, N, Q)
            4. Basic (H, K, R)
            5. Hydrophobic (I, L, M, V)
            6. Aromatic (F, W, Y)
        Hydrophobic-polar maps to a mere two categories:
            1. Hydrophobic (A, F, G, I, L, M, P, V, W, Y)
            2. Polar (C, D, E, H, K, N, Q, R, S, T)
    verbose : bool
        Whether or not to print a lot of stuff
    jaccard_threshold : float
        Value between 0 and 1. By default, the (empirically-chosen) "best"
        threshold is chosen for each alphabet. For "protein" and  "dayhoff",
        the default is 0.5, and for "hydrophobic-polar," it is 0.8, since it is
        so lossy it's more likely to match random sequence. These thresholds
        were determined empirically with a pre-chosen human RNA-seq dataset and
        human peptides.
    description : str
        The identifier in the sequence file, i.e. the name or descriptor of the
        sequence
    noncoding_file_handle : None or file
        If not None, write noncoding nucleotide reads to this file handle
    coding_nucleotide_file_handle : None or file
        If not None, write coding nucleotides reads to this file handle
    low_complexity_peptide_file_handle : None or file
        If not None, write low complexity peptide sequences to this file handle

    Returns
    -------
    max_fraction_in_peptide_db : float
        Of all reading frames, the maximum number of k-mers that matches the
        peptide database
    max_n_kmers: int
        Of all reading frames, the maximum number of k-mers observed in the
        translated, encoded peptide
    special_case : str or None
        Additional message to write in the output csv describing the reason
        why this sequence is or isn't protein-coding
    """
    # Convert to BioPython sequence object for translation
    seq = Seq(sequence)

    # In case this is used from the Python API and the default threshold isn't
    # specified
    jaccard_threshold = get_jaccard_threshold(jaccard_threshold, alphabet)

    # Convert to BioPython sequence object for translation
    translations = six_frame_translation_no_stops(seq)
    # For all translations, use the one with the maximum number of k-mers
    # in the databse
    max_n_kmers = 0
    max_fraction_in_peptide_db = 0
    if len(translations) == 0:
        return SingleReadScore(
            np.nan, np.nan, PROTEIN_CODING_CATEGORIES['stop_codons'])

    translations = {
        frame: translation
        for frame, translation in translations.items()
        if len(translation) > peptide_ksize
    }
    if len(translations) == 0:
        return SingleReadScore(
            np.nan,
            np.nan,
            PROTEIN_CODING_CATEGORIES['too_short_peptide'])

    for frame, translation in translations.items():
        # Convert back to string
        translation = str(translation)

        # Maybe reencode to dayhoff/hp space
        encoded = encode_peptide(translation, alphabet)

        is_kmer_low_complexity, n_kmers = evaluate_is_kmer_low_complexity(
            encoded, peptide_ksize)

        if is_kmer_low_complexity:
            maybe_write_fasta(description + f" translation_frame: {frame}",
                              low_complexity_peptide_file_handle, translation)
            category = LOW_COMPLEXITY_CATEGORIES[alphabet]
            return SingleReadScore(np.nan, n_kmers, category)

        fraction_in_peptide_db, n_kmers = score_single_translation(
            encoded,
            peptide_bloom_filter,
            peptide_ksize,
            molecule=alphabet,
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
            seqname = f'{description} translation_frame: {frame} ' \
                      f'jaccard: {fraction_in_peptide_db}'
            write_fasta(sys.stdout, seqname, translation)
            maybe_write_fasta(seqname, coding_nucleotide_file_handle, sequence)

    if max_fraction_in_peptide_db <= jaccard_threshold:
        maybe_write_fasta(description, noncoding_file_handle, sequence)
    return SingleReadScore(max_fraction_in_peptide_db, max_n_kmers, None)


def maybe_write_fasta(description, file_handle, sequence):
    """Write fasta to file handle if it is not None"""
    if file_handle is not None:
        write_fasta(file_handle, description, sequence)


def score_reads(reads,
                peptide_bloom_filter,
                jaccard_threshold=None,
                molecule='protein',
                verbose=False,
                coding_nucleotide_fasta=None,
                noncoding_nucleotide_fasta=None,
                low_complexity_nucleotide_fasta=None,
                low_complexity_peptide_fasta=None):
    """Assign a coding score to each read. Where the magic happens."""
    jaccard_threshold = get_jaccard_threshold(jaccard_threshold, molecule)
    peptide_ksize = peptide_bloom_filter.ksize()

    scoring_lines = []
    nucleotide_ksize = 3 * peptide_ksize

    fastas, file_handles = maybe_open_fastas(coding_nucleotide_fasta,
                                             low_complexity_nucleotide_fasta,
                                             low_complexity_peptide_fasta,
                                             noncoding_nucleotide_fasta)
    with screed.open(reads) as records:
        for record in tqdm(records):
            description = record['name']
            sequence = record['sequence']
            if verbose:
                print(description)

            jaccard, n_kmers, special_case = maybe_score_single_read(
                description, fastas, file_handles, jaccard_threshold, molecule,
                nucleotide_ksize, peptide_bloom_filter, peptide_ksize,
                sequence, verbose)

            line = get_coding_score_line(description, jaccard,
                                         jaccard_threshold, n_kmers,
                                         special_case)
            scoring_lines.append(line)

    maybe_close_files(file_handles)

    # Concatenate all the lines into a single dataframe
    scoring_df = pd.DataFrame(scoring_lines, columns=SCORING_DF_COLUMNS)

    # Add the reads that were used to generate these scores as a column
    scoring_df['filename'] = reads
    return scoring_df


def get_jaccard_threshold(jaccard_threshold, molecule):
    if jaccard_threshold is None:
        if molecule == 'hp' or molecule == 'hydrophobic-polar':
            jaccard_threshold = DEFAULT_HP_JACCARD_THRESHOLD
        else:
            jaccard_threshold = DEFAULT_JACCARD_THRESHOLD
    return jaccard_threshold


def maybe_score_single_read(description, fastas, file_handles,
                            jaccard_threshold, molecule, nucleotide_ksize,
                            peptide_bloom_filter, peptide_ksize, sequence,
                            verbose):
    """Check if read is low complexity/too short, otherwise score it"""
    # Check if nucleotide sequence is low complexity
    is_fastp_low_complexity = evaluate_is_fastp_low_complexity(sequence)
    if is_fastp_low_complexity:
        n_kmers = np.nan
        jaccard, n_kmers, special_case = too_short_or_low_complexity_nucleotide(
            description, fastas, n_kmers, sequence)
    else:
        jaccard, n_kmers, special_case = score_single_read(
            sequence,
            peptide_bloom_filter,
            peptide_ksize,
            molecule,
            verbose,
            jaccard_threshold=jaccard_threshold,
            description=description,
            noncoding_file_handle=file_handles['noncoding_nucleotide'],
            coding_nucleotide_file_handle=file_handles['coding_nucleotide'],
            low_complexity_peptide_file_handle=file_handles[
                'low_complexity_peptide'])

        if verbose > 1:
            click.echo(f"Jaccard: {jaccard}, n_kmers = {n_kmers}", err=True)
    return SingleReadScore(jaccard, n_kmers, special_case)


def too_short_or_low_complexity_nucleotide(
        description, fastas, n_kmers, sequence):
    if n_kmers > 0:
        jaccard = np.nan
        special_case = "Low complexity nucleotide"
        maybe_write_fasta(description, fastas['low_complexity_nucleotide'],
                          sequence)
    else:
        jaccard = np.nan
        n_kmers = np.nan
        special_case = 'Read length was shorter than 3 * peptide ' \
                       'k-mer size'
    return SingleReadScore(jaccard, n_kmers, special_case)


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
    fastas = {
        "noncoding_nucleotide": noncoding_nucleotide_fasta,
        "coding_nucleotide": coding_nucleotide_fasta,
        "low_complexity_nucleotide": low_complexity_nucleotide_fasta,
        "low_complexity_peptide": low_complexity_peptide_fasta
    }
    file_handles = {}
    for seqtype, fasta in fastas.items():
        if fasta is not None:
            file_handles[seqtype] = open_and_announce(fasta, seqtype)
        else:
            file_handles[seqtype] = None
    return OutputFastasHandles(fastas, file_handles)


def maybe_write_csv(coding_scores, csv):
    if csv:
        click.echo(f"Writing coding scores of reads to {csv}", err=True)
        coding_scores.to_csv(csv, index=False)


def maybe_write_json_summary(coding_scores, json_summary, filenames,
                             bloom_filter, molecule, peptide_ksize,
                             jaccard_threshold, groupby='filename'):
    if not json_summary:
        # Early exit if json_summary is not True
        return
    empty_coding_categories = make_empty_coding_categories(molecule)

    if coding_scores.empty:
        summary = {
            'input_files': filenames,
            'jaccard_info': {
                "count": 0,
                "mean": None,
                "std": None,
                "min": None,
                "25%": None,
                "50%": None,
                "75%": None,
                "max": None
            },
            'classification_value_counts': empty_coding_categories,
            'classification_percentages': empty_coding_categories,
            'histogram_n_coding_frames_per_read': {
                "1": 0,
                "2": 0,
                "3": 0,
                "4": 0,
                "5": 0,
                "6": 0,
            },
            'histogram_n_coding_frames_per_read_percentages': {
                "1": 0,
                "2": 0,
                "3": 0,
                "4": 0,
                "5": 0,
                "6": 0,
            },
        }
    else:
        summary = generate_coding_summary(coding_scores, bloom_filter,
                                          molecule, peptide_ksize,
                                          jaccard_threshold)
    with open(json_summary, 'w') as f:
        click.echo(f"Writing extract_coding summary to {json_summary}",
                   err=True)
        json.dump(summary, fp=f)
    return summary


def generate_coding_summary(coding_scores, bloom_filter_filename, molecule,
                            peptide_ksize, jaccard_threshold):
    translation_frame_percentages, translation_frame_counts = \
        get_n_translated_frames_per_read(
            coding_scores)

    files = coding_scores.filename.unique().tolist()

    classification_percentages, classification_value_counts = \
        get_n_per_coding_classification(coding_scores, molecule)

    # Get Jaccard distributions, count, min, max, mean, stddev, median
    jaccard_info = coding_scores.jaccard_in_peptide_db.describe() \
        .to_dict()
    summary = assemble_summary_dict(
        bloom_filter_filename, classification_percentages,
        classification_value_counts, translation_frame_counts,
        translation_frame_percentages, files, jaccard_info,
        jaccard_threshold, molecule, peptide_ksize)
    return summary


def make_empty_coding_categories(molecule):
    coding_categories = dict.fromkeys(PROTEIN_CODING_CATEGORIES.values(), 0)
    molecule_low_complexity_key = LOW_COMPLEXITY_CATEGORIES[molecule]
    coding_categories[molecule_low_complexity_key] = 0
    return coding_categories


def assemble_summary_dict(bloom_filter, classification_percentages,
                          classification_value_counts,
                          coding_per_read_histogram_for_json,
                          coding_per_read_histogram_percentages_for_json,
                          files, jaccard_info, jaccard_threshold, molecule,
                          peptide_ksize):
    summary = {
        'input_files': files,
        'jaccard_info': jaccard_info,
        'classification_value_counts':
            classification_value_counts,
        'classification_percentages':
            classification_percentages,
        'histogram_n_coding_frames_per_read':
            coding_per_read_histogram_for_json,
        'histogram_n_coding_frames_per_read_percentages':
            coding_per_read_histogram_percentages_for_json,
        'peptide_bloom_filter': bloom_filter,
        'peptide_alphabet': molecule,
        'peptide_ksize': peptide_ksize,
        'jaccard_threshold': jaccard_threshold,
    }
    return summary


def get_n_per_coding_classification(coding_scores, molecule):
    # Initialize to all zeros
    counts = make_empty_coding_categories(molecule)
    # Replace with observed sequences
    counts.update(
        coding_scores.classification.value_counts().to_dict())
    # Convert to series to make percentage calculation easy
    classifications = pd.Series(counts)

    # Initialize to all zeros
    percentages = make_empty_coding_categories(molecule)
    percentages_series = 100 * classifications / classifications.sum()
    # Replace with observations
    percentages.update(percentages_series.to_dict())
    return PercentagesCounts(percentages, counts)


def get_n_translated_frames_per_read(coding_scores, col='n_translated_frames'):
    """Of all coding sequences, get number of possible translations"""
    predicted_coding = coding_scores.query('classification == "Coding"')

    n_coding_per_read = predicted_coding.read_id.value_counts()
    n_coding_per_read.index = n_coding_per_read.index.astype(str)

    n_coding_per_read.name = col
    n_coding_per_read_df = n_coding_per_read.to_frame()

    # Number of reading frames per read, per filename
    coding_per_read_histogram = n_coding_per_read_df[col].value_counts()
    index = coding_per_read_histogram.index.astype(str)
    coding_per_read_histogram.index = 'Number of reads with ' + index + \
                                      " putative protein-coding translations"

    # Total number of coding reads
    total = coding_per_read_histogram.sum()

    coding_per_read_histogram_percentages = \
        100 * coding_per_read_histogram / total
    histogram_for_json = \
        coding_per_read_histogram.to_dict()
    percentages_for_json = \
        coding_per_read_histogram_percentages.to_dict()
    return PercentagesCounts(percentages_for_json, histogram_for_json)


@click.command()
@click.argument('peptides', nargs=1)
@click.argument('reads', nargs=-1)
@click.option('--peptide-ksize',
              default=None, type=int,
              help="K-mer size of the peptide sequence to use. Defaults for"
                   " different alphabets are, "
                   f"protein: {DEFAULT_PROTEIN_KSIZE}"
                   f", dayhoff: {DEFAULT_DAYHOFF_KSIZE},"
                   f" hydrophobic-polar: {DEFAULT_HP_KSIZE}")
@click.option("--save-peptide-bloom-filter",
              is_flag=True,
              default=False,
              help="If specified, save the peptide bloom filter. "
                   "Default filename is the name of the fasta file plus a "
                   "suffix denoting the protein encoding and peptide ksize")
@click.option('--peptides-are-bloom-filter',
              is_flag=True,
              default=False,
              help="Peptide file is already a bloom filter")
@click.option('--jaccard-threshold',
              default=None, type=click.FLOAT, callback=validate_jaccard,
              help="Minimum fraction of peptide k-mers from read in the "
                   "peptide database for this read to be called a "
                   f"'coding read'. Default: {DEFAULT_JACCARD_THRESHOLD} for"
                   f" protein and dayhoff encodings, and "
                   f"{DEFAULT_HP_JACCARD_THRESHOLD} for hydrophobic-polar "
                   f"(hp) encoding")
@click.option('--alphabet', '--encoding', '--molecule',
              default='protein',
              help="The type of amino acid encoding to use. Default is "
                   "'protein', but 'dayhoff' or 'hydrophobic-polar' can be "
                   "used")
@click.option('--csv',
              default=False,
              help='Name of csv file to write with all sequence reads and '
                   'their coding scores')
@click.option('--json-summary',
              default=False,
              help='Name of json file to write summarization of coding/'
                   'noncoding/other categorizations, the '
                   'min/max/mean/median/stddev of Jaccard scores, and other')
@click.option("--coding-nucleotide-fasta",
              help="If specified, save the coding nucleotides to this file")
@click.option("--noncoding-nucleotide-fasta",
              help="If specified, save the noncoding nucleotides to this file")
@click.option("--low-complexity-nucleotide-fasta",
              help="If specified, save the low-complexity nucleotides to this"
                   " file")
@click.option("--low-complexity-peptide-fasta",
              help="If specified, save the low-complexity peptides to this "
                   "file")
@click.option('--tablesize', type=BASED_INT,
              default="1e8",
              help='Size of the bloom filter table to use')
@click.option('--n-tables', type=int,
              default=DEFAULT_N_TABLES,
              help='Size of the bloom filter table to use')
@click.option("--long-reads",
              is_flag=True,
              help="If set, then only considers reading frames starting with "
                   "start codon (ATG) and ending in a stop codon "
                   "(TAG, TAA, TGA)")
@click.option("--verbose", is_flag=True, help="Print more output")
def cli(peptides,
        reads,
        peptide_ksize=None,
        save_peptide_bloom_filter=True,
        peptides_are_bloom_filter=False,
        jaccard_threshold=None,
        alphabet='protein',
        csv=False,
        json_summary=False,
        coding_nucleotide_fasta=None,
        noncoding_nucleotide_fasta=None,
        low_complexity_nucleotide_fasta=None,
        low_complexity_peptide_fasta=None,
        tablesize=DEFAULT_MAX_TABLESIZE, n_tables=DEFAULT_N_TABLES,
        long_reads=False,
        verbose=False):
    """Writes coding peptides from reads to standard output

    \b
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
    save_peptide_bloom_filter : str or bool
        Whether or not to save the created bloom filter to file. If a string,
        save to this filename
    peptides_are_bloom_filter : bool
        Input ilfe of peptides is already a bloom filter
    jaccard_threshold : float
        Value between 0 and 1. By default, the (empirically-chosen) "best"
        threshold is chosen for each alphabet. For "protein" and  "dayhoff",
        the default is 0.5, and for "hydrophobic-polar," it is 0.8, since it is
        so lossy it's more likely to match random sequence. These thresholds
        were determined empirically with a pre-chosen human RNA-seq dataset and
        human peptides.
    alphabet : str
        One of "protein"|"peptide", "dayhoff", or "hydrophobic-polar"|"hp" to
        encode the protein-coding space. Where "protein"|"peptide" is the
        original 20-letter amino acid encoding, Dayhoff ("dayhoff") is a lossy
        6-letter encoding that categorizes the amino acids into:
            1. Cysteine,
            2. Small (A, G, P, S, T)
            3. Acid and Amide (D, E, N, Q)
            4. Basic (H, K, R)
            5. Hydrophobic (I, L, M, V)
            6. Aromatic (F, W, Y)
        Hydrophobic-polar maps to a mere two categories:
            1. Hydrophobic (A, F, G, I, L, M, P, V, W, Y)
            2. Polar (C, D, E, H, K, N, Q, R, S, T)
    csv : str
        Save the coding scores as a csv to this file
    long_reads : bool -- NOT IMPLEMENTED!!
        Input sequencing reads are long reads. Not implemented, but the plan
        is, instead of doing 6-frame translation as on the short reads, test
        all ATG (start codon) to stop codon reading frames for the one(s) that
        matches the known peptide database best. Unknown whether this requires
        new thresholds
    coding_nucleotide_fasta : None or str
        If specified, save coding nucleotide sequence to this file
    noncoding_nucleotide_fasta : None or str
        If specified, save noncoding nucleotide sequence to this file
    low_complexity_nucleotide_fasta : None or str
        If specified, save low complexity nucleotide sequence to this file
    low_complexity_peptide_fasta : None or str
        If specified, save low complexity peptide sequence to this file
    verbose : bool
        Whether or not to print lots of stuff. Can specify multiple, e.g. -vv
        if you really like having everything in stdout

    \b
    Returns
    -------
    coding_peptides : str
        Outputs a fasta-formatted sequence of translated peptides
    """
    # \b above prevents re-wrapping of paragraphs

    if long_reads:
        raise NotImplementedError("Not implemented! ... yet :)")

    peptide_bloom_filter = maybe_make_peptide_bloom_filter(
        peptides, peptide_ksize, alphabet, peptides_are_bloom_filter,
        n_tables=n_tables, tablesize=tablesize)
    click.echo("\tDone!", err=True)

    if not peptides_are_bloom_filter:
        peptide_bloom_filter_filename = maybe_save_peptide_bloom_filter(
            peptides, peptide_bloom_filter, alphabet,
            save_peptide_bloom_filter)
    else:
        peptide_bloom_filter_filename = peptides

    dfs = []
    for reads_file in reads:
        df = score_reads(
            reads_file,
            peptide_bloom_filter,
            jaccard_threshold=jaccard_threshold,
            molecule=alphabet,
            verbose=verbose,
            coding_nucleotide_fasta=coding_nucleotide_fasta,
            noncoding_nucleotide_fasta=noncoding_nucleotide_fasta,
            low_complexity_nucleotide_fasta=low_complexity_nucleotide_fasta,
            low_complexity_peptide_fasta=low_complexity_peptide_fasta)
        dfs.append(df)

    coding_scores = pd.concat(dfs, ignore_index=True)

    maybe_write_csv(coding_scores, csv)
    maybe_write_json_summary(coding_scores, json_summary, reads,
                             bloom_filter=peptide_bloom_filter_filename,
                             molecule=alphabet, peptide_ksize=peptide_ksize,
                             jaccard_threshold=jaccard_threshold)


if __name__ == '__main__':
    cli()
