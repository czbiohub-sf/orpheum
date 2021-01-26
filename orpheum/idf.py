from collections import Counter
import math


def filter_idf(hashes, idf, mean_idf_per_cell):
    return set(x for x in hashes if idf[x] > mean_idf_per_cell)


def _get_single_cell_mean_idf(signature, idf):
    min_hashes = signature.minhash.get_mins()

    vsum = sum(idf[x] for x in min_hashes)
    mean_idf = vsum / len(min_hashes)
    return mean_idf


def get_mean_idf_per_cell(siglist, idf):
    idf_sums = sum(_get_single_cell_mean_idf(s, idf) for s in siglist)
    mean_idf = idf_sums / len(siglist)

    return mean_idf


def get_term_frequency(hash_abundance):
    """Convert raw term abundance to relative frequency in document"""
    total_counts = sum(hash_abundance.values())
    return {k: v / total_counts for k, v in hash_abundance.items()}


def get_document_frequency(siglist):
    document_frequency = Counter()
    for s in siglist:
        document_frequency.update(s.minhash.get_mins())
    return document_frequency


def get_inverse_document_frequency(siglist):

    # total number of documents
    N = len(siglist)
    document_frequency = get_document_frequency(siglist)
    inverse_document_frequency = {
        k: math.log(N / v) for k, v in document_frequency.items()
    }
    return inverse_document_frequency
