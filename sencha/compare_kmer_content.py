from functools import partial
import itertools
import multiprocessing
import os
from pprint import pprint
import random
from typing import Sequence
import time

import click
import pandas as pd
import screed
from sourmash.logging import notify

# Divergence time estimates in millions of years
# from http://www.timetree.org/ on 2019-08-26
from sencha.sequence_encodings import (
    amino_keto_ize,
    weak_strong_ize,
    purine_pyrimidize,
    encode_peptide,
)

MOLECULES_TO_COMPARE = (
    "peptide20",
    "hsdm17",
    "sdm12",
    "aa9",
    "botvinnik8",
    "dayhoff6",
    "gbmr4",
    "hp2",
)

divergence_estimates = pd.Series(
    {
        "Amniota": 312,
        "Bilateria": 824,
        "Boreoeutheria": 96,
        # Old world monkeys
        "Catarrhini": 29.4,
        "Euarchontoglires": 76,
        # Bony vertebrates
        "Euteleostomi": 435,
        "Eutheria": 105,
        # Jawed vertebrates
        "Gnathostomata": 473,
        # A primate suborder
        "Haplorrhini": 67,
        # Great apes (includes orangutan)
        "Hominidae": 15.8,
        # Gorilla, human, chimp
        "Homininae": 9.1,
        # Apes (includes gibbons)
        "Hominoidea": 20.2,
        "Mammalia": 177,
        "Opisthokonta": 1105,
        "Primates": 74,
        # tetrapods and the lobe-finned fishes
        "Sarcopterygii": 413,
        "Simiiformes": 43,
        # Tetrapods - 4-limbed
        "Tetrapoda": 352,
        # Includes Eutheria (placental mammals) and
        # Metatheria (maruspials)
        "Theria": 159,
        "NA": 0,
    }
)
divergence_estimates = divergence_estimates.sort_values()

KSIZES = (
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    23,
    24,
    25,
)
COLUMNS = "id1", "id2", "ksize", "jaccard"


# Hydrophobic/hydrophilic mapping

# From: Brüne, D., Andrade-Navarro, M. A., & Mier, P. (2018).
# Proteome-wide comparison between the amino acid composition of domains and
# linkers. BMC Research Notes, 1–6. http://doi.org/10.1186/s13104-018-3221-0


def sanitize_id(value):
    """Takes first non-whitespace as ID, replaces pipes with underscore

    Cribbed from https://stackoverflow.com/a/295466/1628971
    """
    value = value.split()[0].replace("|", "__")
    return value


def kmerize(seq, ksize):
    """Return the set of unique k-mers from the sequence"""
    return set(seq[i : i + ksize] for i in range(len(seq) - ksize + 1))


def jaccardize(set1, set2):
    """Compute jaccard index of two sets"""
    denominator = min(len(set1), len(set2))
    if denominator > 0:
        return len(set1.intersection(set2)) / denominator
    else:
        return denominator


def kmerize_and_jaccard(seq1, seq2, ksize, debug=False):
    kmers1 = set(seq1[i : i + ksize] for i in range(len(seq1) - ksize + 1))
    kmers2 = set(seq2[i : i + ksize] for i in range(len(seq2) - ksize + 1))
    jaccard = jaccardize(kmers1, kmers2)
    if debug:
        print("len(kmers1):", len(kmers1))
        print("len(kmers2):", len(kmers2))
        print(f"jaccard: {jaccard}")
    return jaccard


def kmer_comparison_table(id1, seq1, id2, seq2, molecule_name, ksizes=KSIZES):
    lines = []
    for ksize in ksizes:
        jaccard = kmerize_and_jaccard(seq1, seq2, ksize)
        if jaccard > 0:
            line = [id1, id2, ksize, jaccard]
            lines.append(line)
        else:
            # If jaccard=0 at a small ksize, then all future jaccards will also
            # be 0 --> break and exit
            remaining_lines = [[id1, id2, k, 0] for k in range(ksize, max(ksizes) + 1)]
            lines.extend(remaining_lines)
            break
    df = pd.DataFrame(lines, columns=COLUMNS)
    df["alphabet"] = molecule_name
    return df


def compare_peptide_seqs(
    id1_seq1, id2_seq2, ksizes=KSIZES, alphabets=MOLECULES_TO_COMPARE
):
    # Unpack the tuples
    id1, seq1 = id1_seq1
    id2, seq2 = id2_seq2

    dfs = []
    for alphabet in alphabets:
        reencoded1 = encode_peptide(seq1, alphabet)
        reencoded2 = encode_peptide(seq2, alphabet)

        df = kmer_comparison_table(
            id1, reencoded1, id2, reencoded2, molecule_name=alphabet, ksizes=ksizes
        )
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)
    return df


def compare_nucleotide_seqs(id1_seq1, id2_seq2, ksizes=KSIZES):
    # Unpack the tuples
    id1, seq1 = id1_seq1
    id2, seq2 = id2_seq2

    purine_pyrimidine1 = purine_pyrimidize(seq1)
    purine_pyrimidine2 = purine_pyrimidize(seq2)

    purine_primimdine_df = kmer_comparison_table(
        id1,
        purine_pyrimidine1,
        id2,
        purine_pyrimidine2,
        molecule_name="purine_pyrimidine",
        ksizes=ksizes,
    )

    weak_strong1 = weak_strong_ize(seq1)
    weak_strong2 = weak_strong_ize(seq2)

    weak_strong_df = kmer_comparison_table(
        id1, weak_strong1, id2, weak_strong2, molecule_name="weak_strong", ksizes=ksizes
    )

    amino_keto1 = amino_keto_ize(seq1)
    amino_keto2 = amino_keto_ize(seq2)

    amino_keto_df = kmer_comparison_table(
        id1, amino_keto1, id2, amino_keto2, molecule_name="amino_keto", ksizes=ksizes
    )

    nucleotide_df = kmer_comparison_table(
        id1, seq1, id2, seq2, molecule_name="nucleotide", ksizes=ksizes
    )

    df = pd.concat(
        [purine_primimdine_df, nucleotide_df, weak_strong_df, amino_keto_df],
        ignore_index=True,
    )
    return df


def compare_seqs(id1_seq1, id2_seq2, ksizes=KSIZES, moltype="protein"):
    if moltype == "protein":
        return compare_peptide_seqs(id1_seq1, id2_seq2, ksizes)
    elif moltype.lower() == "dna":
        return compare_nucleotide_seqs(id1_seq1, id2_seq2, ksizes)
    else:
        raise ValueError(
            f"{moltype} is not a valid molecule type! Only "
            f"'protein' and 'dna' are supported"
        )


def compare_args_unpack(args, ksizes, moltype):
    """Helper function to unpack the arguments. Written to use in pool.imap as
    it can only be given one argument."""
    return compare_seqs(*args, ksizes=ksizes, moltype=moltype)


def get_comparison_at_index(
    index,
    seqlist1,
    seqlist2=None,
    ksizes=KSIZES,
    n_background=100,
    moltype="protein",
    verbose=False,
    paired_seqlists=True,
    intermediate_csv=False,
    intermediate_parquet=False,
    no_final_concatenation=False,
):
    """Returns similarities of all combinations of seqlist1 seqlist2 at index

    Parameters
    ----------
    index : int
        generate masks from this image
    seqlist1 : list
        List of (id, seq) tuples
    seqlist2 : list, optional (default None)
        List of (id, seq) tuples. If None, then an all-by-all comparison of
        sequences in seqlist1 is performed, as if seqlist1 was provided as
        seqlist2.
    ksizes : iterable of int
        K-mer sizes to extract and compare the sequences on
    moltype : str, optional (default "protein")
        One of "protein" or "dna" -- for knowing which alphabets to use
    verbose : boolean, default False
    n_background : int, optional (default 100)
        When paired_seqlist is True, how many random background sequences to
        choose from seqlist2
    paired_seqlists : bool, optional (default True)
        If True, then seqlist1 and seqlist2 have sequences at the same index
        that need to be compared, i.e. index 0 across the two. Best used when
        seqlist1 and seqlist2 are lists of homologous protein sequences across
        two different species
    intermediate_parquet : bool
        Write intermediate file of all comparisons at index i to an
        IO-efficient parquet format
    intermediate_csv : bool
        Write intermediate file of all comparisons at index i to an
        csv format

    Returns
    -------
    comparison_df_list : list
        list of pandas.DataFrame tables for the combinations of seqlist1 at
        index, compared to seqlist2
    """
    startt = time.time()
    id1 = seqlist1[index][0]
    id1_sanitized = sanitize_id(id1)
    csv = id1_sanitized + ".csv"
    parquet = id1_sanitized + ".parquet"
    if os.path.exists(parquet):
        notify(f"Found {parquet} already exists for {id1}, skipping", end="\r")
        return []
    if os.path.exists(csv):
        notify(f"Found {csv} already exists for {id1}, skipping", end="\r")
        return []

    if seqlist2 is not None:
        if paired_seqlists:
            seq_iterator = get_paired_seq_iterator(
                index, n_background, seqlist1, seqlist2, verbose
            )
        else:
            seq_iterator = itertools.product([seqlist1[index]], seqlist2)
    else:
        seq_iterator = itertools.product([seqlist1[index]], seqlist1[index + 1 :])

    func = partial(compare_args_unpack, ksizes=ksizes, moltype=moltype)
    comparision_df_list = list(map(func, seq_iterator))
    notify(
        "comparison for index {} (id: {}) done in {:.5f} seconds",
        index,
        id1,
        time.time() - startt,
        end="\n",
    )

    if intermediate_csv or intermediate_parquet:
        df = pd.concat(comparision_df_list)
        if intermediate_csv:
            df.to_csv(csv)
        if intermediate_parquet:
            df.to_parquet(parquet)
        del df
    if no_final_concatenation:
        del comparision_df_list
        return []
    else:
        return comparision_df_list


def get_paired_seq_iterator(index, n_background, seqlist1, seqlist2, verbose):
    """Compare index i for seqlist1 and seqlist2, plus background

    Choose random seqs from seqlist2 for background

    Parameters
    ----------
    index : int
        Position in `seqlist1` and `seqlist2` to compare to one another, e.g.
        when positions in `seqlist1` and `seqlist2` are known homologs of one
        another
    n_background: int
        Since the pairs are assumed to be matched, we need a background of the
        overall seq-seq similarity, so this assigns the number of random pairs
        chosen in `seqlist2` to be compared to the sequence at `index`
    seqlist2 : list
        List of (id, seq) tuples
    seqlist1 : list
        List of (id, seq) tuples
    verbose : bool
        If True, then print the number of foreground and background pairs to
        be created by this iterator

    Returns
    -------
    seq_iterator : iterable
        Generator of ((seq1, id1), (seq2, id2)) pairs
    """
    pairs_iterator = [(seqlist1[index], seqlist2[index])]
    random_seqlist2 = random.sample(seqlist2, n_background)
    this_index_seqlist1 = [seqlist1[index]] * n_background
    background_pairs = list(zip(this_index_seqlist1, random_seqlist2))
    if verbose:
        print("background_pairs:")
        pprint(background_pairs)
    seq_iterator = list(itertools.chain(*[pairs_iterator, background_pairs]))
    if verbose:
        print("seq_iterator:")
        pprint(seq_iterator)
    return seq_iterator


def compare_all_seqs(
    seqlist1,
    seqlist2=None,
    n_jobs=4,
    ksizes=KSIZES,
    moltype="protein",
    n_background=100,
    paired_seqlists=True,
    intermediate_csv=False,
    intermediate_parquet=False,
    no_final_concatenation=False,
):
    """Compare k-mer content of sequences across k-mer sizes and alphabets

    Parameters
    ----------
    seqlist1 : list
        List of (id, seq) tuples
    seqlist2 : list, optional
        List of (id, seq) tuples. If None, then an all-by-all comparison of
        sequences in seqlist1 is performed, as if seqlist1 was provided as
        seqlist2.
    ksizes : iterable of int
        K-mer sizes to extract and compare the sequences on
    moltype : str
        One of "protein" or "dna" -- for knowing which alphabets to use
    n_background : int
        When paired_seqlist is True, how many random background sequences to
        choose from seqlist2
    n_jobs : int
        Number of jobs for multiprocessing
    paired_seqlists : bool
        If True, then seqlist1 and seqlist2 have sequences at the same index
        that need to be compared, i.e. index 0 across the two. Best used when
        seqlist1 and seqlist2 are lists of homologous protein sequences across
        two different species
    intermediate_parquet : bool
        Write intermediate file of all comparisons at index i to an
        IO-efficient parquet format
    intermediate_csv : bool
        Write intermediate file of all comparisons at index i to an
        csv format

    Returns
    -------
    kmer_comparisons : pandas.DataFrame
        A table of seq1_id, seq2_id, ksize, alphabet encoding, jaccard
        similarity

    Raises
    ------
    ValueError:
        If paired_seqlist=True and seqlist1 and seqlist2 are of different
        lengths, as the comparison is done pairwise across both, as if the
        'zip' operator was used.

    """
    if seqlist2 is not None:
        if paired_seqlists and len(seqlist1) != len(seqlist2):
            raise ValueError(
                "When comparing pairs of sequences, can only "
                "compare two sequences of equal length"
            )
        elif not paired_seqlists:
            # Want seqlist1 to be shorter so that there are fewer, bigger jobs
            # to minimize thread spawning costs
            if len(seqlist2) > len(seqlist1):
                # Swap the seqlist orders so seqlist1 is the shorter one
                old_seqlist1 = seqlist1
                old_seqlist2 = seqlist2
                seqlist2 = old_seqlist1
                seqlist1 = old_seqlist2
    else:
        seqlist2 = seqlist1

    n = len(seqlist1)
    m = len(seqlist2)

    n_comparisons = n * m
    t0 = time.time()
    len_seqlist1 = len(seqlist1)
    notify(f"Number of comparisons: {n} * {m} = {n_comparisons:,}")

    # Initialize the function using func.partial with the common arguments like
    # siglist, ignore_abundance, downsample, for computing all the signatures
    # The only changing parameter that will be mapped from the pool is the
    # index
    func = partial(
        get_comparison_at_index,
        seqlist1=seqlist1,
        seqlist2=seqlist2,
        n_background=n_background,
        ksizes=ksizes,
        moltype=moltype,
        paired_seqlists=paired_seqlists,
        intermediate_csv=intermediate_csv,
        intermediate_parquet=intermediate_parquet,
        no_final_concatenation=no_final_concatenation,
    )
    notify("Created similarity func")

    # Initialize multiprocess.pool
    pool = multiprocessing.Pool(processes=n_jobs)

    # Calculate chunk size, by default pool.imap chunk size is 1
    chunksize, extra = divmod(len_seqlist1, n_jobs)
    if extra:
        chunksize += 1
    notify("Calculated chunk size for multiprocessing")

    # This will not generate the results yet, since pool.imap returns a
    # generator
    result = pool.imap(func, range(len_seqlist1), chunksize=chunksize)
    notify("Initialized multiprocessing pool.imap")

    peptide_kmer_comparisons = pd.concat(itertools.chain(*result), ignore_index=True)

    notify(f"Total time: {time.time() - t0}")
    return peptide_kmer_comparisons


@click.command()
@click.argument("fastas", nargs=-1)
@click.option(
    "--fastas2",
    default=None,
    help="Optional. Instead of doing an all-by-all comparison of "
    "the provided fasta arguments, do fastas2 vs fastas "
    "arguments",
)
@click.option(
    "--alphabets",
    default=",".join(MOLECULES_TO_COMPARE),
    help="Which protein-coding alphabet to use for comparisons",
)
@click.option(
    "--ksize-min", default=2, type=click.INT, help="k-mer sizes to use for comparison"
)
@click.option(
    "--ksize-max", default=25, type=click.INT, help="k-mer sizes to use for comparison"
)
@click.option(
    "--ksize-step", default=1, type=click.INT, help="k-mer sizes to use for comparison"
)
@click.option(
    "--parquet",
    default=None,
    help="If provided, save table to a space-efficient and fast-IO "
    "Apache Parquet format file of this name. This format is "
    "compatible with Python/Pandas and R.",
)
@click.option(
    "--no-csv", is_flag=True, default=False, help="Don't output csv to stdout"
)
@click.option(
    "--paired-seqlists",
    is_flag=True,
    default=False,
    help="Fastas are paired 1:1 by item, so compare item 1 in "
    "fasta1 to item 1 in fasta2",
)
@click.option(
    "--intermediate-parquet",
    is_flag=True,
    default=False,
    help="If provided, write a parquet file for each individual "
    "comparison at index i, in current directory",
)
@click.option(
    "--intermediate-csv",
    is_flag=True,
    default=False,
    help="If provided, write a csv file for each individual "
    "comparison at index i, in current directory",
)
@click.option(
    "--no-final-concatenation",
    is_flag=True,
    default=False,
    help="If provided, then don't concatenate ",
)
@click.option(
    "--processes",
    "-p",
    default=2,
    type=click.INT,
    help="Number of processes to use for parallelization",
)
def cli(
    fastas,
    fastas2,
    alphabets,
    ksize_min,
    ksize_max,
    ksize_step,
    parquet,
    no_csv=False,
    paired_seqlists=False,
    intermediate_csv=False,
    intermediate_parquet=False,
    no_final_concatenation=False,
    processes=2,
):
    """Compute k-mer similarity of all pairwise sequences"""
    if len(fastas) == 0:
        raise ValueError(
            "No sequence files provided! " "Argument 'fastas' is required!"
        )
    seqlist = parse_fastas(fastas)
    if fastas2 is not None:
        seqlist2 = parse_fastas([fastas2])
    else:
        seqlist2 = None

    if no_final_concatenation:
        if not intermediate_parquet or not intermediate_csv:
            raise Exception(
                "--no-final-concatenation provided but neither "
                "--intermediate-parquet nor --intermediate-csv "
                "provided, so no output will be created!"
            )

    # add 1 to max since range is not inclusive of last interval
    ksizes = list(range(ksize_min, ksize_max + 1, ksize_step))

    comparisons = compare_all_seqs(
        seqlist,
        seqlist2=seqlist2,
        n_jobs=processes,
        ksizes=ksizes,
        moltype="protein",
        paired_seqlists=paired_seqlists,
        intermediate_csv=intermediate_csv,
        intermediate_parquet=intermediate_parquet,
        no_final_concatenation=no_final_concatenation,
    )

    # Only write the final output if there is a final concatenation
    if not no_final_concatenation:
        if parquet is not None:
            comparisons.to_parquet(parquet)
        if not no_csv:
            click.echo(comparisons.to_csv(index=False))


def parse_fastas(fastas: Sequence):
    """Open fasta files and create list of (id, seq) tuples

    Parameters
    ----------
    fastas : Sequence
        List of strings or paths of fastas to open with Screed

    Returns
    -------
    seqlist : List[Tuple[str, str]]
        List of (id, seq) tuple

    """
    seqlist = []
    for fasta in fastas:
        with screed.open(fasta) as records:
            for record in records:
                seq_id = record["name"]
                seq = record["sequence"]
                seqlist.append((seq_id, seq))
    if len(seqlist) == 0:
        raise ValueError(f"No sequences found in files: {' '.join(fastas)}!")
    return seqlist


if __name__ == "__main__":
    cli()
