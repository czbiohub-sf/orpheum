
from collections import Counter, defaultdict
import math


def filter_idf(hashes, idf, mean_idf_per_cell):
    return set(x for x in hashes if idf[x] > mean_idf_per_cell)


def _get_single_cell_mean_idf(signature, idf):
    min_hashes = signature.minhash.get_mins()
    
    vsum = sum(idf[x] for x in min_hashes)
    mean_idf = vsum/len(min_hashes)
    return mean_idf


def get_mean_idf_per_cell(siglist, idf):
    idf_sums = sum(_get_single_cell_mean_idf(s, idf) for s in siglist)
    mean_idf = idf_sums/len(siglist)
    
    return mean_idf


def get_term_frequency(hash_abundance):
    """Convert raw term abundance to relative frequency in document"""
    total_counts = sum(hash_abundance.values())
    return {k: v/total_counts for k, v in hash_abundance.items()}


def get_document_frequency(siglist):
    document_frequency = Counter()
    for s in siglist:
        document_frequency.update(s.minhash.get_mins())
    return document_frequency


def get_inverse_document_frequency(siglist):
    
    # total number of documents
    N = len(siglist)
    document_frequency = get_document_frequency(siglist)
    inverse_document_frequency = {k: math.log(N/v) for k, v in document_frequency.items()}
    return inverse_document_frequency


def make_siglist_groups(siglist, groups):
    """Map groups names in series to lists of signatures

    Returns
    -------
    siglist_grouped : dict
        dictionary of string name to a list of signatures
    """
    groups_siglist = [groups[x.name()] for x in siglist]

    siglist_grouped = defaultdict(list)
    for animal, sig in zip(groups_siglist, siglist):
        siglist_grouped[animal].append(sig)
    return siglist_grouped


def perform_tf_idf(siglist, groups):
    siglist_grouped = make_siglist_groups(siglist, groups)

    siglist_grouped_tf_idf = {}

    siglist_species_tf_idf = []

    for species, species_sigs in siglist_grouped.items():
        species_sigs = list(species_sigs)
        species_idf = get_inverse_document_frequency(species_sigs)
        species_mean_idf_per_cell = idf.get_mean_idf_per_cell(species_sigs,
                                                              species_idf)
        species_sigs_filtered_idf = [idf.filter_idf(x.minhash.get_mins(),
                                                    species_idf,
                                                    species_mean_idf_per_cell)
                                     for x in species_sigs]
        siglist_grouped_tf_idf[species] = species_sigs_filtered_idf
        siglist_species_tf_idf.extend(species_sigs_filtered_idf)
