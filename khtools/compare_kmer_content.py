from functools import partial
import itertools
import multiprocessing
from pprint import pprint
import random
import time


import pandas as pd
from sourmash.logging import notify

# Divergence time estimates in millions of years
# from http://www.timetree.org/ on 2019-08-26
from .sequence_encodings import amino_keto_ize, \
    weak_strong_ize, purine_pyrimidize, encode_peptide

MOLECULES_TO_COMPARE = 'aa20', 'dayhoff6', 'hp2', 'botvinnik8', 'aa9', \
                       'gbmr4', 'sdm12', 'hsdm17'

divergence_estimates = pd.Series({"Amniota": 312,
                                  'Bilateria': 824,
                                  "Boreoeutheria": 96,

                                  # Old world monkeys
                                  'Catarrhini': 29.4,
                                  "Euarchontoglires": 76,

                                  # Bony vertebrates
                                  'Euteleostomi': 435,
                                  'Eutheria': 105,

                                  # Jawed vertebrates
                                  'Gnathostomata': 473,

                                  # A primate suborder
                                  'Haplorrhini': 67,

                                  # Great apes (includes orangutan)
                                  'Hominidae': 15.8,

                                  # Gorilla, human, chimp
                                  'Homininae': 9.1,

                                  # Apes (includes gibbons)
                                  'Hominoidea': 20.2,

                                  'Mammalia': 177,
                                  "Opisthokonta": 1105,
                                  'Primates': 74,

                                  # tetrapods and the lobe-finned fishes
                                  'Sarcopterygii': 413,
                                  'Simiiformes': 43,

                                  # Tetrapods - 4-limbed
                                  'Tetrapoda': 352,

                                  # Includes Eutheria (placental mammals) and
                                  # Metatheria (maruspials)
                                  'Theria': 159,

                                  'NA': 0})
divergence_estimates = divergence_estimates.sort_values()


KSIZES = 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, \
    21, 23, 24, 25
COLUMNS = 'id1', 'id2', 'ksize', 'jaccard'

# Hydrophobic/hydrophilic mapping

# From: Brüne, D., Andrade-Navarro, M. A., & Mier, P. (2018).
# Proteome-wide comparison between the amino acid composition of domains and
# linkers. BMC Research Notes, 1–6. http://doi.org/10.1186/s13104-018-3221-0


def kmerize(seq, ksize):
    """Return the set of unique k-mers from the sequence"""
    return set(seq[i:i + ksize] for i in range(len(seq) - ksize + 1))


def jaccardize(set1, set2):
    """Compute jaccard index of two sets"""
    denominator = min(len(set1), len(set2))
    if denominator > 0:
        return len(set1.intersection(set2)) / denominator
    else:
        return denominator


def kmerize_and_jaccard(seq1, seq2, ksize, debug=False):
    kmers1 = set(seq1[i:i + ksize] for i in range(len(seq1) - ksize + 1))
    kmers2 = set(seq2[i:i + ksize] for i in range(len(seq2) - ksize + 1))
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
        line = [id1, id2, ksize, jaccard]
        lines.append(line)
    df = pd.DataFrame(lines, columns=COLUMNS)
    df['alphabet'] = molecule_name
    return df


def compare_peptide_seqs(id1_seq1, id2_seq2, ksizes=KSIZES,
                         alphabets=MOLECULES_TO_COMPARE):
    # Unpack the tuples
    id1, seq1 = id1_seq1
    id2, seq2 = id2_seq2

    dfs = []
    for alphabet in alphabets:
        reencoded1 = encode_peptide(seq1, alphabet)
        reencoded2 = encode_peptide(seq2, alphabet)

        df = kmer_comparison_table(id1, reencoded1, id2, reencoded2,
                                   molecule_name=alphabet,
                                   ksizes=ksizes)
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
        id1, purine_pyrimidine1, id2, purine_pyrimidine2,
        molecule_name='purine_pyrimidine', ksizes=ksizes)

    weak_strong1 = weak_strong_ize(seq1)
    weak_strong2 = weak_strong_ize(seq2)

    weak_strong_df = kmer_comparison_table(
        id1, weak_strong1, id2, weak_strong2,
        molecule_name='weak_strong', ksizes=ksizes)

    amino_keto1 = amino_keto_ize(seq1)
    amino_keto2 = amino_keto_ize(seq2)

    amino_keto_df = kmer_comparison_table(
        id1, amino_keto1, id2, amino_keto2,
        molecule_name='amino_keto', ksizes=ksizes)

    nucleotide_df = kmer_comparison_table(id1, seq1, id2, seq2,
                                          molecule_name='nucleotide',
                                          ksizes=ksizes)

    df = pd.concat([purine_primimdine_df, nucleotide_df,
                    weak_strong_df, amino_keto_df], ignore_index=True)
    return df


def compare_seqs(id1_seq1, id2_seq2, ksizes=KSIZES, moltype='protein'):
    if moltype == 'protein':
        return compare_peptide_seqs(id1_seq1, id2_seq2, ksizes)
    elif moltype.lower() == 'dna':
        return compare_nucleotide_seqs(id1_seq1, id2_seq2, ksizes)


def compare_args_unpack(args, ksizes, moltype):
    """Helper function to unpack the arguments. Written to use in pool.imap as it
    can only be given one argument."""
    return compare_seqs(*args, ksizes=ksizes, moltype=moltype)


def get_comparison_at_index(index, seqlist1, seqlist2,
                            ksizes=KSIZES, n_background=100,
                            moltype='protein', verbose=False):
    """Returns similarities of all the combinations of signature at index in the
    siglist with the rest of the indices starting at index + 1. Doesn't
    redundantly calculate signatures with all the other indices prior to
    index - 1

    :param int index: generate masks from this image
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate a distance metric
        based on the cosine similarity.
    :param boolean downsample by max_hash if True
    :param siglist list of signatures
    :return: list of similarities for the combinations of signature at index
    with rest of the signatures from index+1
    """
    startt = time.time()
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
    func = partial(compare_args_unpack, ksizes=ksizes, moltype=moltype)
    comparision_df_list = list(map(func, seq_iterator))
    notify(
        "comparison for index {} done in {:.5f} seconds",
        index,
        time.time() - startt,
        end='\r')
    return comparision_df_list


def compare_all_seqs(seqlist1, seqlist2=None, n_jobs=4, ksizes=KSIZES,
                     moltype='protein', n_background=100):
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
    if seqlist2 is not None:
        if len(seqlist1) != len(seqlist2):
            raise ValueError("Can only compare two sequences of equal length")

    t0 = time.time()
    len_seqlist1 = len(seqlist1)

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
        moltype=moltype)
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

    peptide_kmer_comparisons = pd.concat(
        itertools.chain(*result), ignore_index=True)

    notify(f"Total time: {time.time() - t0}")
    return peptide_kmer_comparisons
