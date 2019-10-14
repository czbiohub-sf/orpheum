

def jaccard_sigs(i, j, siglist):
    return siglist[i].jaccard(siglist[j])


def jaccard(sample1, sample2):
    """Jaccard similarity between two sets"""
    intersection = len(sample1.intersection(sample2))
    union = len(sample1.union(sample2))
    return intersection / union
