
import jaccard_utils

import itertools

import numpy as np
from scipy.spatial.distance import squareform


def _compare_serial(siglist, iterator):
    n = len(siglist)
    values = np.ones((n, n))
    
    for i, j in iterator:
        jaccard = jaccard_utils.jaccard_sigs(i, j, siglist)
        
        values[i, j] = jaccard
        values[j, i] = jaccard
        
    return values


def compare_all_pairs(siglist, n_jobs=None):
    n = len(siglist)
    iterator = itertools.combinations(range(n), 2)
    
    
    if n_jobs is None:
        values = _compare_serial(siglist, iterator)
    else:
        # This creates a condensed distance matrix
        condensed = Parallel(n_jobs=n_jobs)(delayed(jaccard_utils.jaccard_sigs)(i, j, siglist) for i, j in iterator)
        values = squareform(condensed)
        
    return values


def get_similarity_difference(similarity, groupby=('ksize', 'molecule', 'num_hashes')):
    """Calculate difference in similarity from "true" aka maximum sampling similarity
    
    Parameters
    ----------
    similiarity : pandas.DataFrame
        tidy dataframe with "cell1", "cell2", "similarity" columns,
        plus groupby columns
    
    """
    CELL_INDEX = ['cell1', 'cell2']

    max_num_hashes = similarity.num_hashes.max()
    true_similarity = similarity.query('num_hashes == @max_num_hashes')
    true_similarity = true_similarity.set_index(CELL_INDEX).sort_index()

    test_similarity = similarity.query('num_hashes != @max_num_hashes')
    test_similarity = test_similarity.set_index(CELL_INDEX).sort_index()

    similarity_difference = test_similarity.groupby(groupby).apply(
        lambda x: x['similarity'].subtract(true_similarity['similarity']))
    similarity_difference = similarity_difference.unstack(level=[0, 1, 2]).reset_index()
    similarity_difference = similarity_difference.rename(columns={0: 'similarity_difference'})
    
    return similarity_difference