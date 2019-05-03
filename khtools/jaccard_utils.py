import os
import tempfile

from joblib import Parallel, delayed, load, dump
from itertools import combinations
from collections import defaultdict



from .idf import filter_idf

def jaccard_sigs(i, j, siglist):
    return siglist[i].jaccard(siglist[j])

def jaccard_sigs_idf(i, j, siglist, idf, mean_idf_per_cell):
    i_hashes = filter_idf(siglist[i].get_mins(), idf, mean_idf_per_cell)
    j_hashes = filter_idf(siglist[j].get_mins(), idf, mean_idf_per_cell)
    return jaccard(i_hashes, j_hashes)


def jaccard(sample1, sample2):
    """Jaccard similarity between two sets"""
    intersection = len(sample1.intersection(sample2))
    union = len(sample1.union(sample2))
    return intersection/union


def estimate_jaccard(a, b, normalize=True):
    """Estimate jaccard similarity between A and B.

    Normalizes A and B to be the same size, and takes the sketch of A U B
    """
    min_length = min(len(a), len(b))

    # Renormalize the sketch lengths to be identical
    if normalize:
        a = set(sorted(a)[:min_length])
        b = set(sorted(b)[:min_length])

    sketch_a_union_b = set(sorted(a.union(b))[:min_length])

    numerator = len(sketch_a_union_b.intersection(a).intersection(b))
    denominator = len(sketch_a_union_b)
    try:
        return numerator / denominator
    except ZeroDivisionError:
        return 0


def memmap_siglist(siglist):
    """Write a memory-mapped array of signatures"""
    temp_folder = tempfile.mkdtemp()
    filename = os.path.join(temp_folder, 'siglist.mmap')
    if os.path.exists(filename): os.unlink(filename)
    _ = dump(siglist, filename)
    large_memmap = load(filename, mmap_mode='r+')
    return large_memmap


def jaccard_sets_parallel(siglist, n_jobs=16):
    """Estimate jaccard similarity between sets of values"""
    memmapped = memmap_siglist(siglist)
    values_idf = Parallel(n_jobs=n_jobs, require='sharedmem',
                          backend='threading')(
        delayed(estimate_jaccard)(x, y) for x, y in combinations(memmapped, 2))
    return values_idf



def jaccard_tf_idf(siglist, series):
    series_groups = [series[x.name()] for x in siglist]

    siglist_grouped = defaultdict(list)
    for animal, sig in zip(series_groups, siglist):
        siglist_grouped[animal].append(sig)
