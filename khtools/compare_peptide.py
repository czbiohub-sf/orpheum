from functools import partial
import itertools
import math
import multiprocessing
import time


import pandas as pd
from sourmash.logging import notify



KSIZES = 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25
COLUMNS = 'id1', 'id2', 'ksize', 'jaccard'

AMINO_ACID_SINGLE_LETTERS = "R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"

DAYHOFF_MAPPING = {
    "C": "a",

    # Small
    "A": "b",
    "G": "b",
    "P": "b",
    "S": "b",
    "T": "b",

    # Acid and amide
    "D": "c",
    "E": "c",
    "N": "c",
    "Q": "c",

    # Basic
    "H": "d",
    "K": "d",
    "R": "d",

    # Hydrophobic
    "I": "e",
    "L": "e",
    "M": "e",
    "V": "e",

    # Aromatic
    "F": "f",
    "W": "f",
    "Y": "f"
}

## Hydrophobic/hydrophilic mapping
HP_MAPPING = {
    # Hydrophobic
    "A": "h",
    "F": "h",
    "G": "h",
    "I": "h",
    "L": "h",
    "M": 'h',
    "P": "h",
    "V": "h",
    "W": "h",
    "Y": "h",

    # Hydrophilic - polar
    "N": 'p',
    "C": 'p',
    "S": "p",
    "T": "p",
    "D": "p",
    "E": "p",
    "R": "p",
    "H": "p",
    "K": "p",
    "Q": "p"
}

# From: Brüne, D., Andrade-Navarro, M. A., & Mier, P. (2018).
# Proteome-wide comparison between the amino acid composition of domains and
# linkers. BMC Research Notes, 1–6. http://doi.org/10.1186/s13104-018-3221-0
BOTVINNIK_MAPPING = {
    # Small and hydrophobic
    "A": "a",
    "G": "a",

    # Hydrophobic
    "L": "b",
    "I": "b",
    "V": "b",

    # Aromatic, not W
    "F": "c",
    "Y": "c",

    # Polar or charged
    # Phosphorylate-able
    "S": "d",
    "T": "d",

    # Polar, uncharged
    "N": "e",
    "Q": "e",

    # Polar, negatively charged
    "D": "f",
    "E": "f",

    # Polar, positively charged
    # Not histidine
    "R": "g",
    "K": "g",

    # Special
    "C": "h",
    "M": "i",
    "W": "j",
    "H": "k",
    "Q": "l",
    "P": "m"
}

assert all(x in DAYHOFF_MAPPING for x in AMINO_ACID_SINGLE_LETTERS)
assert all(x in HP_MAPPING for x in AMINO_ACID_SINGLE_LETTERS)
assert all(x in BOTVINNIK_MAPPING for x in AMINO_ACID_SINGLE_LETTERS)

DAYHOFF_TRANSLATION = str.maketrans(DAYHOFF_MAPPING)

HP_TRANSLATION = str.maketrans(HP_MAPPING)

BOTVINNIK_TRANSLATION = str.maketrans(BOTVINNIK_MAPPING)

def dayhoffize(seq):
    return seq.translate(DAYHOFF_TRANSLATION)


def hpize(seq):
    return seq.translate(HP_TRANSLATION)

def botvinnikize(seq):
    return seq.translate(BOTVINNIK_TRANSLATION)


def kmerize(seq, ksize):
    return set(seq[i:i+ksize] for i in range(len(seq)-ksize+1))

def jaccardize(set1, set2):
    denominator = min(len(set1), len(set2))
    if denominator > 0:
        return len(set1.intersection(set2))/denominator
    else:
        return denominator

def kmerize_and_jaccard(seq1, seq2, ksize, debug=False):
    kmers1 = set(seq1[i:i+ksize] for i in range(len(seq1)-ksize+1))
    kmers2 = set(seq2[i:i+ksize] for i in range(len(seq2)-ksize+1))
    jaccard = jaccardize(kmers1, kmers2)
    if debug:
        print("len(kmers1):", len(kmers1))
        print("len(kmers2):", len(kmers2))
        print(f"jaccard: {jaccard}")
    return jaccard


def kmer_comparison_table(id1, seq1, id2, seq2, molecule, ksizes=KSIZES):
    lines = []
    for ksize in ksizes:
        jaccard = kmerize_and_jaccard(seq1, seq2, ksize)
        line = [id1, id2, ksize, jaccard]
        lines.append(line)
    df = pd.DataFrame(lines, columns=COLUMNS)
    df['molecule'] = molecule
    return df


def nCr(n, r):
    f = math.factorial
    return f(n) // (f(r) * f(n - r))

def compare_peptide_seqs(id1_seq1, id2_seq2, ksizes=KSIZES):
    # Unpack the tuples
    id1, seq1 = id1_seq1
    id2, seq2 = id2_seq2

    protein_df = kmer_comparison_table(id1, seq1, id2, seq2,
                                       molecule='protein', ksizes=ksizes)

    botvinnik1 = botvinnikize(seq1)
    botvinnik2 = botvinnikize(seq2)

    botvinnik_df = kmer_comparison_table(id1, botvinnik1, id2, botvinnik2,
                                  molecule='botvinnik', ksizes=ksizes)

    dayhoff1 = dayhoffize(seq1)
    dayhoff2 = dayhoffize(seq2)

    dayhoff_df = kmer_comparison_table(id1, dayhoff1, id2, dayhoff2,
                                       molecule='dayhoff', ksizes=ksizes)

    hp1 = hpize(seq1)
    hp2 = hpize(seq2)

    hp_df = kmer_comparison_table(id1, hp1, id2, hp2,
                                  molecule='hydrophobic-polar', ksizes=ksizes)

    df = pd.concat([protein_df, botvinnik_df, dayhoff_df, hp_df], ignore_index=True)
    return df


def compare_args_unpack(args, ksizes):
    """Helper function to unpack the arguments. Written to use in pool.imap as it
    can only be given one argument."""
    return compare_peptide_seqs(*args, ksizes=ksizes)


def get_comparison_at_index(index, seqlist, ksizes=KSIZES):
    """Returns similarities of all the combinations of signature at index in the
    siglist with the rest of the indices starting at index + 1. Doesn't redundantly
    calculate signatures with all the other indices prior to index - 1

    :param int index: generate masks from this image
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate a distance metric
        based on the cosine similarity.
    :param boolean downsample by max_hash if True
    :param siglist list of signatures
    :return: list of similarities for the combinations of signature at index with
    rest of the signatures from index+1
    """
    startt = time.time()
    seq_iterator = itertools.product([seqlist[index]], seqlist[index + 1:])
    func = partial(compare_args_unpack, ksizes=ksizes)
    comparision_df_list = list(map(func, seq_iterator))
    notify(
        "comparison for index {} done in {:.5f} seconds",
        index,
        time.time() - startt,
        end='\r')
    return comparision_df_list


def compare_all_seqs(seqlist, n_jobs=4, ksizes=KSIZES):
    """
    
    Parameters
    ----------
    seqlist : list
        List of (id, seq) tuples
    n_jobs : int
        Number of jobs for multiprocessing

    Returns
    -------

    """
    t0 = time.time()
    length_seqlist = len(seqlist)

    # Initialize the function using func.partial with the common arguments like
    # siglist, ignore_abundance, downsample, for computing all the signatures
    # The only changing parameter that will be mapped from the pool is the index
    func = partial(
        get_comparison_at_index,
        seqlist=seqlist,
        ksizes=ksizes)
    notify("Created similarity func")

    # Initialize multiprocess.pool
    pool = multiprocessing.Pool(processes=n_jobs)

    # Calculate chunk size, by default pool.imap chunk size is 1
    chunksize, extra = divmod(length_seqlist, n_jobs)
    if extra:
        chunksize += 1
    notify("Calculated chunk size for multiprocessing")

    # This will not generate the results yet, since pool.imap returns a generator
    result = pool.imap(func, range(length_seqlist), chunksize=chunksize)
    notify("Initialized multiprocessing pool.imap")

    peptide_kmer_comparisons = pd.concat(itertools.chain(*result), ignore_index=True)

    notify(f"Total time: {time.time() - t0}")
    return peptide_kmer_comparisons
