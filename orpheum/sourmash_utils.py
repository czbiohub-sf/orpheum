import itertools
import json
import warnings

from matplotlib.colors import rgb2hex
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
import seaborn as sns


from .extract_metadata import extract_cell_metadata

keys_for_length = ("mins",)

keys_for_values = ("ksize", "alphabet")

keys_to_print = keys_for_length + keys_for_values

BLADDER_CELL_ID = "10X_P4_3_GCGAGAACACATGGGA"
TISSUE_CHANNEL = "bladder-10X_P4_3_"

PALETTE_NAMES = (
    "tab10",
    "Set2",
    "tab20",
    "Set3",
    "tab20b",
    "Accent",
    "tab20c",
    "Dark2",
    "Paired",
    "Pastel1",
    "Set1",
    "Pastel2",
)

METADATA_COL = ["cell_ontology_class", "free_annotation"]


def _describe_single(signature):
    data = []
    for sig in signature["signatures"]:
        d = dict(name=signature["name"])
        for key, value in sig.items():
            if key not in keys_to_print:
                continue

            if key in keys_for_length:
                v = len(value)
                k = "n_" + key
            elif key in keys_for_values:
                v = value
                k = key
            d[k] = v
        data.append(d)
    return data


def describe(filename):
    with open(filename) as f:
        signature = json.load(f)
    data = itertools.chain(*map(_describe_single, signature))

    description = pd.DataFrame(list(data))
    description["log10_n_mins"] = np.log10(description.n_mins)

    return description


def get_single_cell(cell_id, matrix, metadata, name, ksize, ignore_abundance):
    cell = matrix.loc[:, cell_id].to_frame()
    cell.columns = ["similarity"]
    cell.index = metadata.index
    cell = cell.join(metadata)
    cell = cell.sort_values("similarity", ascending=False)

    cell["name"] = name
    cell["ksize"] = ksize
    cell["ignore_abundance"] = ignore_abundance

    return cell


def get_unique_ordered_categories(categories):
    unique = sorted(
        categories.dropna().unique(),
        key=lambda x: x.lower() if isinstance(x, str) else x,
    )
    return unique


def single_category_colors(categories, palette):
    """Convert a series of categories to colors"""
    no_na = categories.dropna()
    item_to_color = pd.Series(index=categories.index)

    unique = get_unique_ordered_categories(no_na)
    n_unique = len(unique)
    colors = [rgb2hex(x) for x in sns.color_palette(palette, n_colors=n_unique)]
    category_to_color = dict(zip(unique, colors))

    data = [category_to_color[c] for c in no_na]
    item_to_color[no_na.index] = data
    return item_to_color


def category_colors(categories, palette):
    """Convert a dataframe of categories to colors

    Params
    ------
    categories : pandas.DataFrame or pandas.Series
    palette : string or dict
        If a string, then "categories" must be a Series and this only labels
        one series
        If a dict, then corresponds to the column names of "categories"
        dataframe
    """
    item_to_color = None
    if isinstance(categories, pd.DataFrame):
        if palette is None:
            palette = dict(zip(categories.columns, PALETTE_NAMES))
        data = [
            single_category_colors(categories[col], palette[col]) for col in categories
        ]
        item_to_color = pd.concat(data, axis=1)
        item_to_color.columns = categories.columns
    elif isinstance(categories, pd.Series):
        item_to_color = single_category_colors(categories, palette)
        item_to_color.name = categories.name

    item_to_color = item_to_color.fillna("#262626")
    return item_to_color


BETWEENS = ("rows", "cols")


def calculate_linkage(data, metric, method, between="cols"):
    assert between in BETWEENS

    if between == "cols":
        data = data.T

    D = pdist(data, metric)
    try:
        from fastcluster import linkage
    except ImportError:
        from scipy.cluster.hierarchy import linkage

        warnings.warn(
            "'fastcluster' not installed! This may take a while ... "
            "Install with 'pip install fastcluster'"
        )

    Z = linkage(D, method)

    try:
        from polo import optimal_leaf_ordering

        optimal_Z = optimal_leaf_ordering(Z, D)
    except ImportError:
        warnings.warn(
            "'polo' not installed! Dendrogram will not be optimal "
            "leaf ordered. Install with 'pip install polo'"
        )
        optimal_Z = Z

    return optimal_Z


def plaidplot(
    data,
    row_categories=None,
    col_categories=None,
    row_palette=None,
    col_palette=None,
    metric="euclidean",
    method="ward",
    xticklabels=[],
    yticklabels=[],
    cmap="GnBu",
    **kwargs,
):

    col_linkage = calculate_linkage(data, metric, method, between="cols")
    row_linkage = calculate_linkage(data, metric, method, between="rows")

    row_colors = (
        category_colors(row_categories, row_palette)
        if "row_colors" not in kwargs
        else kwargs.pop("row_colors")
    )
    col_colors = (
        category_colors(col_categories, col_palette)
        if "col_colors" not in kwargs
        else kwargs.pop("col_colors")
    )

    if "vmax" not in kwargs:
        kwargs["vmax"] = data.replace(1, np.nan).max().max()

    g = sns.clustermap(
        data,
        col_linkage=col_linkage,
        row_linkage=row_linkage,
        row_colors=row_colors,
        col_colors=col_colors,
        cmap=cmap,
        xticklabels=xticklabels,
        yticklabels=yticklabels,
        **kwargs,
    )
    return g


def plaidplot_square(
    data, metadata, metadata_col="cell_ontology_class", palette="tab20", **kwargs
):
    categories = metadata[metadata_col]
    #     palette = category_colors(categories, palette)

    vmax = data.replace(1, np.nan).max().max()

    return plaidplot(
        data,
        col_palette=palette,
        row_palette=palette,
        col_categories=categories,
        row_categories=categories,
        vmax=vmax,
        **kwargs,
    )


def facet_distplot(df, x="similarity", hue="cell_ontology_class", palette="tab20"):
    hue_order = get_unique_ordered_categories(df[hue])
    g = sns.FacetGrid(df, hue=hue, hue_order=hue_order, palette=palette, height=3)

    vmax = df[x].replace(1, np.nan).max().max()
    # add a bit of left and right padding so it looks nice
    buffer = vmax * 0.05

    g.map(sns.distplot, x, hist=False, kde_kws=dict(shade=True))
    g.add_legend()
    g.set(xlim=(-buffer, vmax + buffer))
    return g


def plaidplot_and_distplot(
    data,
    metadata,
    name,
    ksize,
    ignore_abundance,
    molecule,
    metadata_col=METADATA_COL,
    cell_id=BLADDER_CELL_ID,
    tissue_channel=TISSUE_CHANNEL,
    palette="tab20",
    **kwargs,
):
    plaidplot_grid = plaidplot_square(
        data, metadata, metadata_col=metadata_col, **kwargs
    )
    fig_prefix = (
        f"{tissue_channel}_{name}_k{ksize}_{molecule}_"
        "ignore-abundance={ignore_abundance}"
    )
    png = f"../figures/{fig_prefix}_clustermap.png"
    plaidplot_grid.ax_col_dendrogram.set(title=fig_prefix)
    plaidplot_grid.savefig(png, dpi=300)

    if "row_palette" in kwargs:
        if "cell_ontology_class" in kwargs["row_palette"]:
            palette = kwargs["row_palette"]["cell_ontology_class"]

    df = get_single_cell(cell_id, data, metadata, ksize, name, ignore_abundance)

    distplot_grid = facet_distplot(df, palette=palette, hue=metadata_col)
    pdf = f"../figures/{fig_prefix}_cell={cell_id}_distplot.pdf"
    distplot_grid.ax.set(title=fig_prefix)
    distplot_grid.savefig(pdf)
    return df


def _assemble_metadata(compare, cell_ids, metadata_cols):
    colon_separated = compare.columns.str.contains(":")

    if colon_separated.sum() > 0:
        metadata_colon_separated = extract_cell_metadata(
            compare.columns[colon_separated],
            pattern="(?P<column>\\w+):(?P<value>[\\w-]+)",
        )
        metadata_colon_separated = metadata_colon_separated.drop(
            columns=["cell"], errors="ignore"
        )
        metadata_colon_separated.index = cell_ids[colon_separated]
    else:
        metadata_colon_separated = pd.DataFrame()

    if len(colon_separated) > 0 and not colon_separated.all():
        # metadata is pipe-separated
        # 'leukocyte|Lung|3-F-56|10X_P7_8_GAACATCTCTTGAGGT'

        # input data has to be a list of lists
        data = list(compare.columns[~colon_separated].str.split("|").values)
        metadata_not_colon_separated = pd.DataFrame(data, columns=metadata_cols)
        metadata_not_colon_separated = metadata_not_colon_separated.dropna(
            subset=["cell_id"]
        )
        metadata_not_colon_separated = metadata_not_colon_separated.set_index("cell_id")
    else:
        metadata_not_colon_separated = pd.DataFrame()

    metadata = pd.concat(
        [metadata_not_colon_separated, metadata_colon_separated],
        sort=False,
        ignore_index=False,
    )
    metadata = metadata.drop(columns=["cell_id"], errors="ignore")
    metadata = metadata.reindex(index=cell_ids)
    metadata["method"] = metadata.index.map(lambda x: "10x" if "10X" in x else "FACS")
    # Replace underscores with spaces for cell ontology names for consistency
    metadata["cell_ontology_class"] = metadata["cell_ontology_class"].str.replace(
        "_", " "
    )
    return metadata


def read_compare(
    csv,
    pattern="(?P<column>\\w+):(?P<value>[\\w-]+)",
    metadata_cols=["cell_ontology_class", "tissue", "mouse_id", "cell_id"],
):

    compare = pd.read_csv(csv)

    cell_ids = compare.columns.str.split("|").str[-1]
    cell_ids = cell_ids.str.split(":").str[-1]
    cell_ids = pd.Index(cell_ids, name="cell_id")

    metadata = _assemble_metadata(compare, cell_ids, metadata_cols)

    # Rename to valid cell ids
    compare.index = compare.columns = cell_ids

    compare = compare.reindex(index=metadata.index, columns=metadata.index)

    return compare, metadata


def filter_siglist(siglist, ksize, moltype):
    if moltype == "protein":

        def molfilter(x):
            return x.minhash.is_protein

    else:

        def molfilter(x):
            return not x.minhash.is_protein

    return [s for s in siglist if molfilter(s) and (s.minhash.ksize == ksize)]
