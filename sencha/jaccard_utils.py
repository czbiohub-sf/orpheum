from .idf import filter_idf


def jaccard_sigs(i, j, siglist):
    return siglist[i].jaccard(siglist[j])


def jaccard_sigs_idf(
    i,
    j,
    siglist,
    idf,
    mean_idf_per_cell,
):
    i_hashes = filter_idf(siglist[i].get_mins(), idf, mean_idf_per_cell)
    j_hashes = filter_idf(siglist[j].get_mins(), idf, mean_idf_per_cell)
    return jaccard(i_hashes, j_hashes)


def jaccard(sample1, sample2):
    """Jaccard similarity between two sets"""
    intersection = len(sample1.intersection(sample2))
    union = len(sample1.union(sample2))
    return intersection / union
