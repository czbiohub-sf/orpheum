from itertools import combinations
from collections import Counter, defaultdict
import math


import pandas as pd

from .jaccard_utils import jaccard_sigs_parallel


COMPARISON_COLS = 'animal', 'tissue', 'replicate'
GROUPBY_COLS = [f'same_{x}' for x in COMPARISON_COLS]

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


def get_tf_idf(siglist):
    inverse_document_frequencies = get_inverse_document_frequency(siglist)

    # TODO: Use x.minhash.get_mins(with_abundance=True) for track_abundance=True
    term_frequencies = [1/len(x.minhash.get_mins()) for x in siglist]



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

def tidify_values_idf(values_idf, index, idf_quantile):
    tidy = pd.Series(values_idf, index=index)
    tidy = tidy.reset_index()
    tidy = tidy.rename(
        columns={'level_0': 'sample_1', 'level_1': 'sample_2',
                 0: 'similarity'})
    tidy['idf_quantile'] = idf_quantile
    return tidy


def parallel_many_tf_idf(siglist, groups, idf_quantiles, metadata,
                         comparison_cols=COMPARISON_COLS,
                         n_jobs=16, plot=True):
    siglist_grouped = make_siglist_groups(siglist, groups)
    species_idfs = {species: get_inverse_document_frequency(list(sigs))
                   for species, sigs in siglist_grouped.items()}

    names = [x.name() for x in siglist]
    index = pd.MultiIndex.from_tuples(combinations(names, 2))

    dfs = []
    for idf_quantile in idf_quantiles:
        print(f'idf_quantile: {idf_quantile}')
        values_idf = tf_idf_and_jaccard(siglist_grouped, species_idfs,
                                        idf_quantile, n_jobs)
        tidy = tidify_values_idf(values_idf, index, idf_quantile)
        tidy = add_comparison_cols(tidy, metadata, comparison_cols)
        if plot:
            title = f"IDF Quantile: {idf_quantile}"
            plot_comparison_similarities(tidy, comparison_cols,
                                         title=title)
        dfs.append(tidy)
    similarities_idf_tidy = pd.concat(dfs, ignore_index=True)
    return similarities_idf_tidy




def tf_idf_and_jaccard(siglist_grouped, species_idfs, idf_quantile=0.5,
                       n_jobs=16):

    siglist_species_tf_idf = []

    for species, species_sigs in siglist_grouped.items():
        species_idf = species_idfs[species]
        species_minimum_idf = pd.Series(species_idf).quantile(idf_quantile)
        species_sigs_filtered_idf = [filter_idf(x.minhash.get_mins(),
                                                    species_idf,
                                                    species_minimum_idf)
                                     for x in species_sigs]
        siglist_species_tf_idf.extend(species_sigs_filtered_idf)
    values_idf = jaccard_sigs_parallel(siglist_species_tf_idf, n_jobs=n_jobs)
    return values_idf



def add_comparison_cols(tidy, metadata, comparison_cols=COMPARISON_COLS):
    """

    Parameters
    ----------
    tidy : pandas.DataFrame
        Dataframe of similarities between sample_1, sample_2
    metadata : pandas.DataFrame
        Per-sample metadata for samples in tidy

    Returns
    -------

    """
    tidy_metadata = tidy.join(metadata, on='sample_1') \
        .join(metadata, on='sample_2', rsuffix='_sample_2')

    for col in comparison_cols:
        tidy_metadata[f'same_{col}'] = \
            tidy_metadata[col] == tidy_metadata[f'{col}_sample_2']

        return tidy_metadata


def plot_comparison_similarities(similarity_tidy_metadata,
                                 comparison_cols=COMPARISON_COLS,
                                 title=None):

    for col in comparison_cols:
        groupby_col = 'same_' + col
        g = sns.FacetGrid(similarity_tidy_metadata, hue=groupby_col)
        g.map(sns.distplot, 'similarity')
        g.add_legend()
        g.set(title=None)
