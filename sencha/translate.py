"""
translate.py

Partition reads into coding, noncoding, and low-complexity bins
"""
import csv
import os
import sys
from functools import partial

import click
import numpy as np
import screed
import multiprocessing
from sourmash._minhash import hash_murmur
from sencha.log_utils import get_logger
from sencha.sequence_encodings import encode_peptide
from sencha.compare_kmer_content import kmerize
from sencha.create_save_summary import CreateSaveSummary
from sencha.index import (
    maybe_make_peptide_bloom_filter,
    maybe_save_peptide_bloom_filter,
)
import sencha.constants_index as constants_index
import sencha.fasta_utils as fasta_utils
import sencha.constants_translate as constants_translate
from sencha.translate_single_seq import TranslateSingleSeq


logger = get_logger(__file__)


def get_jaccard_threshold(jaccard_threshold, alphabet):
    if jaccard_threshold is None:
        if alphabet == "hp" or alphabet == "hydrophobic-polar":
            jaccard_threshold = constants_translate.DEFAULT_HP_JACCARD_THRESHOLD
        else:
            jaccard_threshold = constants_translate.DEFAULT_JACCARD_THRESHOLD
    return jaccard_threshold


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
        raise click.BadParameter(
            "--jaccard-threshold needs to be a number"
            " between 0 and 1, but {} was provided".format(value)
        )


def evaluate_is_fastp_low_complexity(sequence):
    """Use fastp's definition of complexity
    By this definition, low complexity sequence is defined by consecutive runs
    of same base in a row, e.g.
    CCCCCCCCCACCACCACCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCACACACCCCCAACAC
    is low complexity. The threshold is 0.3 as used in the fastp prpject:
    https://github.com/OpenGene/fastp
    Parameters
    ----------
    sequence : str
        Sequence to compute complexity on
    complexity_threshold : float, defaault 0.3
        Value between 0 and 1. The default is 0.3 because that is the default
        in the command line program fastp
    Returns
    -------
    is_low_complexity : bool
        Whether or not the sequence passes the complexity threshold
    """
    assert len(sequence) != 0, "sequence is {}".format(sequence)
    complexity = compute_fastp_complexity(sequence)
    return complexity < constants_translate.COMPLEXITY_THRESHOLD


def compute_fastp_complexity(sequence):
    assert len(sequence) != 0, "sequence is {}".format(sequence)
    n_different_consecutively = sum(
        1 for i in range(len(sequence) - 1) if sequence[i] != sequence[i + 1]
    )
    complexity = n_different_consecutively / len(sequence)
    return complexity


def evaluate_is_kmer_low_complexity(sequence, ksize):
    """Check if sequence is low complexity, i.e. mostly repetitive
    By this definition, the sequence is not complex if its number of unique
    k-mers is smaller than half the number of expected k-mers
    """
    try:
        kmers = kmerize(sequence, ksize)
    except ValueError:
        # k-mer size is larger than sequence
        return None, None
    n_kmers = len(kmers)
    complexity = compute_kmer_complexity(sequence, ksize)
    is_low_complexity = n_kmers <= complexity
    return is_low_complexity


def compute_kmer_complexity(sequence, ksize):
    n_possible_kmers_on_sequence = len(sequence) - ksize + 1
    complexity = n_possible_kmers_on_sequence / 2
    return complexity


def score_single_translation(translation, alphabet, peptide_ksize, verbose=False):
    """Score a single translation based on
    fraction of kmers in peptide bloom filter"""
    encoded = encode_peptide(translation, alphabet)
    kmers = list(set(kmerize(str(encoded), peptide_ksize)))
    hashes = [hash_murmur(kmer) for kmer in kmers]
    n_kmers = len(kmers)
    n_kmers_in_peptide_db = sum(1 for h in hashes if peptide_bloom_filter.get(h) > 0)
    if verbose:
        logger.info("\ttranslation: \t".format(encoded))
        logger.info("\tkmers:", " ".join(kmers))

    if verbose:
        kmers_in_peptide_db = {
            (k, h): peptide_bloom_filter.get(h) for k, h in zip(kmers, hashes)
        }
        # Print keys (kmers) only
        logger.info("\tK-mers in peptide database:")
        logger.info(kmers_in_peptide_db)

    fraction_in_peptide_db = n_kmers_in_peptide_db / n_kmers

    return fraction_in_peptide_db, n_kmers


def check_peptide_content(
    description,
    sequence,
    fastas,
    alphabet,
    peptide_ksize,
    jaccard_threshold,
    verbose=False,
):
    """Predict whether a nucleotide sequence could be protein-coding"""

    translations = TranslateSingleSeq(sequence, verbose).six_frame_translation()

    for frame, translation in translations.items():
        if "*" in translation:
            yield constants_translate.SingleReadScore(
                np.nan,
                np.nan,
                constants_translate.PROTEIN_CODING_CATEGORIES["stop_codons"],
                frame,
            )
        elif len(translation) <= peptide_ksize:
            yield constants_translate.SingleReadScore(
                np.nan,
                np.nan,
                constants_translate.PROTEIN_CODING_CATEGORIES["too_short_peptide"],
                frame,
            )
        else:
            encoded = encode_peptide(str(translation), alphabet)
            fraction_in_peptide_db, n_kmers = score_single_translation(
                encoded, alphabet, peptide_ksize, verbose
            )
            is_kmer_low_complex = evaluate_is_kmer_low_complexity(
                encoded, peptide_ksize
            )
            if is_kmer_low_complex:
                fasta_utils.maybe_write_fasta(
                    fastas["low_complexity_peptide"],
                    description + " translation_frame: {}.format(frame)",
                    translation,
                )
                yield constants_translate.SingleReadScore(
                    np.nan,
                    n_kmers,
                    constants_translate.LOW_COMPLEXITY_CATEGORIES[alphabet],
                    frame,
                )
            else:
                if verbose:
                    logger.info(
                        "\t{} is above {}".format(translation, jaccard_threshold)
                    )
                sequencename = "{} translation_frame: {} ".format(
                    description, frame
                ) + "jaccard: {}".format(fraction_in_peptide_db)
                if fraction_in_peptide_db > jaccard_threshold:
                    fasta_utils.maybe_write_fasta(
                        fastas["coding_peptide"], sequencename, translation
                    )
                    fasta_utils.maybe_write_fasta(
                        fastas["coding_nucleotide"], sequencename, sequence
                    )
                    yield constants_translate.SingleReadScore(
                        fraction_in_peptide_db,
                        n_kmers,
                        constants_translate.PROTEIN_CODING_CATEGORIES["coding"],
                        frame,
                    )
                else:
                    fasta_utils.maybe_write_fasta(
                        fastas["noncoding_nucleotide"], sequencename, sequence
                    )
                    yield constants_translate.SingleReadScore(
                        fraction_in_peptide_db,
                        n_kmers,
                        constants_translate.PROTEIN_CODING_CATEGORIES["non_coding"],
                        frame,
                    )


def check_nucleotide_content(description, n_kmers, sequence, fastas):
    """If passes, then this read can move on to checking protein translations

    Evaluates if this reads' nucleotide content doesn't
    pass thresholds to be
    checked for protein-coding-ness
    """
    if n_kmers > 0:
        jaccard = np.nan
        special_case = "Low complexity nucleotide"
        fasta_utils.maybe_write_fasta(
            fastas["low_complexity_nucleotide"], description, sequence
        )
    else:
        jaccard = np.nan
        n_kmers = np.nan
        special_case = "Read length was shorter than 3 * peptide " "k-mer size"
        for i in [1, 2, 3, -1, -2, -3]:
            yield constants_translate.SingleReadScore(jaccard, n_kmers, special_case, i)


def get_coding_score_line(
    description, jaccard, n_kmers, special_case, frame, jaccard_threshold
):
    if special_case is not None:
        line = [description, jaccard, n_kmers, special_case, frame]
    elif jaccard > jaccard_threshold:
        line = [description, jaccard, n_kmers, "Coding", frame]
    else:
        line = [description, jaccard, n_kmers, "Non-coding", frame]
    return line


def maybe_score_single_read(
    alphabet, peptide_ksize, jaccard_threshold, current_reads_file, verbose, fasta_path
):
    """Check if read is low complexity/too short, otherwise score it"""
    fasta_prefix = fasta_path.replace(".fasta", "")
    coding_peptide_fasta = fasta_prefix + "_coding_peptide.fasta"
    csv_name = coding_peptide_fasta.replace(".fasta", ".csv")
    fastas = {
        "noncoding_nucleotide": fasta_prefix + "_" + "noncoding_nucleotide.fasta",
        "coding_nucleotide": fasta_prefix + "_" + "coding_nucleotide.fasta",
        "low_complexity_nucleotide": fasta_prefix
        + "_"
        + "low_complexity_nucleotide.fasta",
        "low_complexity_peptide": fasta_prefix + "_" + "low_complexity_peptide.fasta",
        "coding_peptide": coding_peptide_fasta,
    }
    fastas = fasta_utils.maybe_open_fastas(fastas)
    lines = []
    # writing the fields
    with screed.open(fasta_path) as records:
        for record in records:
            description = record["name"]
            sequence = record["sequence"]
            assert len(sequence) != 0, "sequence is {}".format(sequence)
            # Check if nucleotide sequence is low complexity
            if evaluate_is_fastp_low_complexity(sequence):
                for (
                    jaccard,
                    n_kmers,
                    special_case,
                    frame,
                ) in check_nucleotide_content(description, np.nan, sequence, fastas):
                    if verbose:
                        logger.info(
                            "Jaccard: {}, n_kmers = {} for frame {}".format(
                                jaccard, n_kmers, frame
                            )
                        )

                    line = get_coding_score_line(
                        description,
                        jaccard,
                        n_kmers,
                        special_case,
                        frame,
                        jaccard_threshold,
                    )
                    line.append(current_reads_file)
                    lines.append(line)
            else:
                for (jaccard, n_kmers, special_case, frame,) in check_peptide_content(
                    description,
                    sequence,
                    fastas,
                    alphabet,
                    peptide_ksize,
                    jaccard_threshold,
                    verbose,
                ):
                    if verbose:
                        logger.info(
                            "Jaccard: {}, n_kmers = {} for frame {}".format(
                                jaccard, n_kmers, frame
                            )
                        )

                    line = get_coding_score_line(
                        description,
                        jaccard,
                        n_kmers,
                        special_case,
                        frame,
                        jaccard_threshold,
                    )
                    line.append(current_reads_file)
                    lines.append(line)
    fasta_utils.maybe_close_fastas(fastas)
    # writing the data rows
    with open(csv_name, "w", newline="") as csvfile:
        # creating a csv writer object
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(lines)


class Translate:
    def __init__(self, args):
        """Constructor"""
        self.args = args
        for key in args:
            setattr(self, key, args[key])
        if self.long_reads:
            raise NotImplementedError("Not implemented! ... yet :)")
        self.jaccard_threshold = get_jaccard_threshold(
            self.jaccard_threshold, self.alphabet
        )
        # Use global variable so the functions that are inside multiprocessing
        # don't try to access the nodegraph by self.peptide_bloom_filter
        # and hence avoiding the pickling error while serialization of the NodeGraph by multiprocessing
        # which doesn't exist for nodegraph currently
        # Con is that global variable definition inside a class defintion
        # is advised to be avoided and is said that it might be enforced as incorrect
        # NB PV: This is been said in the docs since Python2.7, but it hasn't changed, so I
        # suggest keeping the global variable declaration here for now
        # https://docs.python.org/3.8/reference/simple_stmts.html#the-global-statement
        # 1) One solution is that if we can call create and declare peptide_bloom_filter inside the
        # if__name__=main function but that would mean switching to argparse from click
        # 2) Another solution is to make nodegraph serializable
        global peptide_bloom_filter
        peptide_bloom_filter = maybe_make_peptide_bloom_filter(
            self.peptides,
            self.peptide_ksize,
            self.alphabet,
            self.peptides_are_bloom_filter,
            n_tables=self.n_tables,
            tablesize=self.tablesize,
        )
        self.verbose = True
        logger.info("\tDone making peptide_bloom_filter!")

        if not self.peptides_are_bloom_filter:
            self.peptide_bloom_filter_filename = maybe_save_peptide_bloom_filter(
                self.peptides,
                peptide_bloom_filter,
                self.alphabet,
                self.save_peptide_bloom_filter,
            )
        else:
            self.peptide_bloom_filter_filename = self.peptides
        self.peptide_ksize = peptide_bloom_filter.ksize()
        self.nucleotide_ksize = 3 * self.peptide_ksize
        self.coding_scores = []
        if not os.path.exists(os.path.abspath(self.intermediate_directory)):
            os.makedirs(self.intermediate_directory)

    def get_jaccard_threshold(self):
        return get_jaccard_threshold(self.jaccard_threshold, self.alphabet)

    def score_reads_per_file(self, reads):
        """Assign a coding score to each read. Where the magic happens."""
        num_records = 0
        with screed.open(reads) as records:
            for record in records:
                num_records += 1
        if num_records == 0:
            return []
        chunksize = fasta_utils.calculate_chunksize(num_records, self.processes)
        fasta_files_split = fasta_utils.split_fasta_files(
            reads, chunksize, self.intermediate_directory
        )
        pool = multiprocessing.Pool(processes=self.processes)
        func = partial(
            maybe_score_single_read,
            self.alphabet,
            self.peptide_ksize,
            self.get_jaccard_threshold(),
            reads,
            self.verbose,
        )
        pool.map(func, fasta_files_split)
        pool.close()
        pool.join()

        # Combine from each of the split fasta's 4 different fastas above to one
        # as in the above file handles
        for key, file_handle in self.file_handles.items():
            for split in fasta_files_split:
                fasta_prefix = split.replace(".fasta", "")
                fasta_path = fasta_prefix + "_" + key + ".fasta"
                with screed.open(fasta_path) as records:
                    for record in records:
                        description = record["name"]
                        sequence = record["sequence"]
                        fasta_utils.maybe_write_fasta(
                            file_handle, description, sequence
                        )

        # Combine and print coding_peptide fastas
        for split in fasta_files_split:
            fasta_prefix = split.replace(".fasta", "")
            fasta_path = fasta_prefix + "_" + "coding_peptide.fasta"
            with screed.open(fasta_path) as records:
                for record in records:
                    description = record["name"]
                    sequence = record["sequence"]
                    fasta_utils.write_fasta(sys.stdout, description, sequence)

        # combine and return scores
        for split in fasta_files_split:
            csv_file = split.replace(".fasta", "_coding_peptide.csv")
            with open(csv_file) as csvfile:
                read_csv = csv.reader(csvfile, delimiter=",")
                for row in read_csv:
                    assert len(row) == 6, "row is {}".format(row)
                    if len(row) == 6:
                        try:
                            line = [
                                str(row[0]),
                                float(row[1]),
                                float(row[2]),
                                str(row[3]),
                                int(row[4]),
                                str(row[5]),
                            ]
                            self.coding_scores.append(line)
                        except:
                            assert 1 == 0, "problems in casting ros {}".format(row)

    def set_coding_scores_all_files(self):
        fastas = {
            "noncoding_nucleotide": self.noncoding_nucleotide_fasta,
            "coding_nucleotide": self.coding_nucleotide_fasta,
            "low_complexity_nucleotide": self.low_complexity_nucleotide_fasta,
            "low_complexity_peptide": self.low_complexity_peptide_fasta,
        }
        self.file_handles = fasta_utils.maybe_open_fastas(fastas)
        for reads_file in self.reads:
            self.score_reads_per_file(reads_file)
        fasta_utils.maybe_close_fastas(self.file_handles)

    def get_coding_scores_all_files(self):
        return self.coding_scores


@click.command()
@click.argument("peptides", nargs=1)
@click.argument("reads", nargs=-1)
@click.option(
    "--peptide-ksize",
    default=None,
    type=int,
    help="K-mer size of the peptide sequence to use. Defaults for"
    " different alphabets are, "
    "protein: {}, dayhoff {}, hydrophobic-polar {}".format(
        constants_index.DEFAULT_PROTEIN_KSIZE,
        constants_index.DEFAULT_DAYHOFF_KSIZE,
        constants_index.DEFAULT_HP_KSIZE,
    ),
)  # noqa
@click.option(
    "--save-peptide-bloom-filter",
    is_flag=True,
    default=False,
    help="If specified, save the peptide bloom filter. "
    "Default filename is the name of the fasta file plus a "
    "suffix denoting the protein encoding and peptide ksize",
)
@click.option(
    "--peptides-are-bloom-filter",
    is_flag=True,
    default=False,
    help="Peptide file is already a bloom filter",
)
@click.option(
    "--jaccard-threshold",
    default=None,
    type=click.FLOAT,
    callback=validate_jaccard,
    help="Minimum fraction of peptide k-mers from read in the "
    "peptide database for this read to be called a "
    "'coding read'.'"
    "Default:"
    "{}".format(constants_translate.DEFAULT_JACCARD_THRESHOLD)
    + " for protein and dayhoff encodings, and "  # noqa
    + "{}".format(constants_translate.DEFAULT_HP_JACCARD_THRESHOLD)  # noqa
    + "for hydrophobic-polar (hp) encoding",  # noqa
)
@click.option(
    "--alphabet",
    "--encoding",
    "--molecule",
    default="protein",
    help="The type of amino acid encoding to use. Default is "
    "'protein', but 'dayhoff' or 'hydrophobic-polar' can be "
    "used",
)
@click.option(
    "--csv",
    default=False,
    help="Name of csv file to write with all sequence reads and " "their coding scores",
)
@click.option(
    "--parquet",
    default=False,
    help="Name of parquet file to write with all sequence reads and "
    "their coding scores",
)
@click.option(
    "--json-summary",
    default=False,
    help="Name of json file to write summarization of coding/"
    "noncoding/other categorizations, the "
    "min/max/mean/median/stddev of Jaccard scores, and other",
)
@click.option(
    "--coding-nucleotide-fasta",
    help="If specified, save the coding nucleotides to this file",
)
@click.option(
    "--noncoding-nucleotide-fasta",
    help="If specified, save the noncoding nucleotides to this file",
)
@click.option(
    "--low-complexity-nucleotide-fasta",
    help="If specified, save the low-complexity nucleotides to this" " file",
)
@click.option(
    "--low-complexity-peptide-fasta",
    help="If specified, save the low-complexity peptides to this " "file",
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
    "--long-reads",
    is_flag=True,
    help="If set, then only considers reading frames starting with "
    "start codon (ATG) and ending in a stop codon "
    "(TAG, TAA, TGA)",
)
@click.option(
    "--intermediate-directory",
    default="/tmp",
    help="If specified this directory is expected to be already created and intermediate files are written here",
)
@click.option(
    "--processes",
    default=1,
    help="number of cores to parallely process on",
)
@click.option("--verbose", is_flag=True, help="Print more output")
def cli(
    peptides,
    reads,
    peptide_ksize=None,
    save_peptide_bloom_filter=True,
    peptides_are_bloom_filter=False,
    jaccard_threshold=None,
    alphabet="protein",
    csv=False,
    parquet=False,
    json_summary=False,
    coding_nucleotide_fasta=None,
    noncoding_nucleotide_fasta=None,
    low_complexity_nucleotide_fasta=None,
    low_complexity_peptide_fasta=None,
    tablesize=constants_index.DEFAULT_MAX_TABLESIZE,
    n_tables=constants_index.DEFAULT_N_TABLES,
    long_reads=False,
    verbose=False,
    processes=1,
    intermediate_directory="/tmp",
):
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
        Input file of peptides is already a bloom filter
    jaccard_threshold : float
        Value between 0 and 1. By default, the (empirically-chosen) "best"
        threshold is chosen for each alphabet. For "protein" and  "dayhoff",
        the default is 0.5, and for "hydrophobic-polar," it is 0.8, since it is
        so lossy it's more likely to match random sequence. These thresholds
        were determined empirically with a pre-chosen human RNA-sequence dataset and
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
    parquet : str
        Save the coding scores as a parquet to this file
    long_reads : bool -- NOT IMPLEMENTED!!
        Input sequenceuencing reads are long reads. Not implemented, but the plan
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
    processes: int
        number of processes to parallel process on
    intermediate_directory: str
        Directory to save intermediate fastas and csv file

    \b
    Returns
    -------
    coding_peptides : str
        Outputs a fasta-formatted sequence of translated peptides
    """
    # \b above prevents re-wrapping of paragraphs
    translate_obj = Translate(locals())
    translate_obj.set_coding_scores_all_files()
    assemble_summary_obj = CreateSaveSummary(
        reads,
        csv,
        parquet,
        json_summary,
        translate_obj.peptide_bloom_filter_filename,
        alphabet,
        translate_obj.peptide_ksize,
        translate_obj.jaccard_threshold,
        translate_obj.coding_scores,
    )
    del translate_obj.coding_scores
    assemble_summary_obj.maybe_write_csv()
    assemble_summary_obj.maybe_write_parquet()
    assemble_summary_obj.maybe_write_json_summary()


if __name__ == "__main__":
    cli()
