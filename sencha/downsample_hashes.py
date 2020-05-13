import time

import pandas as pd

from .sourmash_compare_utils import compare_all_pairs
from .sourmash_utils import filter_siglist


num_hashess = 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000


def downsample_siglist(siglist, downsample_num_hashes=None, downsample_scaled=None):
    if downsample_num_hashes is not None:
        siglist_downsampled = [
            s.minhash.downsample_n(downsample_num_hashes) for s in siglist
        ]
    elif downsample_scaled is not None:
        siglist_downsampled = [
            s.minhash.downsample_scaled(downsample_scaled) for s in siglist
        ]
    else:
        raise ValueError(
            "Either downsample_num_hashes or downsample_scaled must be " "specified!"
        )
    return siglist_downsampled


def compare_downsampled(siglist, num_hashes=None, scaled=None, names=None):
    siglist_downsampled = downsample_siglist(
        siglist, downsample_num_hashes=num_hashes, downsample_scaled=scaled
    )

    if names is None:
        names = [s.name() for s in siglist]
    assert len(names) == len(siglist)

    t0 = time.time()
    values = compare_all_pairs(siglist_downsampled)
    time_delta = time.time() - t0

    df = pd.DataFrame(values, index=names, columns=names).unstack().reset_index()
    df = df.rename(columns={"level_0": "cell1", "level_1": "cell2", 0: "similarity"})

    # Remove all self-similarity
    df = df.query("similarity < 1")
    df.loc[:, "num_hashes"] = num_hashes
    df.loc[:, "time_seconds"] = time_delta
    return df


def subset_ksize_moltype_and_compare_numhashes(
    siglist, ksize, molecule, num_hashes_to_downsample=num_hashess
):
    siglist_filtered = filter_siglist(siglist, ksize, molecule)
    names = [s.name() for s in siglist_filtered]

    dfs = []

    for num_hashes in num_hashes_to_downsample:
        print(f"\tnum_hashes: {num_hashes}")
        df = compare_downsampled(siglist_filtered, num_hashes=num_hashes, names=names)
        dfs.append(df)
    similarity_filtered = pd.concat(dfs)
    print(similarity_filtered.shape)

    similarity_filtered.loc[:, "ksize"] = ksize
    similarity_filtered.loc[:, "alphabet"] = molecule
    return similarity_filtered
