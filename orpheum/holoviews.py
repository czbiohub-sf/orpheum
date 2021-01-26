import itertools
import warnings

# 3rd-party modules
import holoviews as hv
from holoviews import opts, dim
from holoviews.operation.datashader import bundle_graph
import networkx as nx
import pandas as pd

# My handwritten modules
from .s3_utils import savefig
from . import knn
from . import sourmash_utils


# don't warn me about too many figures open
import matplotlib.pyplot as plt

plt.rcParams.update({"figure.max_open_warning": 0})


KSIZES = 9, 12, 15, 21
LOG2SKETCHSIZES = 10, 12, 14, 16
MOLECULES = "dna", "protein"

COLOR_COLS = [
    "species",
    "cell_label",
]

PALETTES = dict(species="Set2", cell_label="tab20")

SKETCH_ID_TEMPLATE = (
    "alphabet-{alphabet}_ksize-{ksize}_" "log2sketchsize-{log2sketchsize}"
)

N_NEIGHBORS = 5


def build_graph_and_plot(
    data,
    metadata,
    n_neighbors,
    color_cols,
    palettes,
    figure_folder,
    figure_prefix,
    title,
):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        graph = knn.nearest_neighbor_graph(
            data,
            metadata,
            n_neighbors=n_neighbors,
            color_cols=color_cols,
            palettes=palettes,
        )

    pos = nx.spring_layout(graph, seed=0)

    for label in color_cols:
        fig, ax = plt.subplots()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            knn.draw_graph(graph, edge_color="black", label_col=label, pos=pos)
        ax.set_title(title)
        figure_suffix = f"graph_nneighbors-{n_neighbors}_colorby-{label}"
        png = f"{figure_folder}/{figure_prefix}_{figure_suffix}.png"
        savefig(fig, png, dpi=150)
    return graph, pos


def get_similarity_graphs(
    csv_template,
    metadata,
    figure_folder,
    groupby="species",
    ksizes=KSIZES,
    log2sketchsizes=LOG2SKETCHSIZES,
    molecules=MOLECULES,
    sketch_id_template=SKETCH_ID_TEMPLATE,
    n_neighbors=N_NEIGHBORS,
    plaidplot=False,
    palettes=PALETTES,
    color_cols=COLOR_COLS,
    verbose=False,
    make_within_groupby_graphs=False,
):
    """Read similarity csvs and create holoviews graphs

    Parameters
    ----------
    csv_template : str
        format-string to insert alphabet, ksize, and log2sketchsize values
        into to get csv. e.g.:
        'similarities_molecule-{alphabet}_ksize-{ksize}_log2sketchsize-{log2sketchsize}.csv'
    metadata : pandas.DataFrame
        Sample-by-feature metadata encoding additional information about
        samples, such as species, cell type label, or tissue
    groupby : str
        Which column of the metadata to groupby to get sub-graphs for
    ksizes : tuple of int
        Which k-mer sizes to look for similarity files for,
        default (9, 12, 15, 21)
    log2sketchsizes : tuple of int
        Which log2 sketch sizes to look for similarity files for,
        default (10, 12, 14, 16)
    molecules : tuple of str
        Which alphabets to use, default both 'dna' and 'protein'
    sketch_id_template : str
        String to use as a unique identifier for the sketch, e.g.
        'alphabet-{alphabet}_ksize-{ksize}_log2sketchsize-{log2sketchsize}'
    plaidplot : bool
        If true, make a clustered heatmap with the sides labeled with the
        color_cols
    palettes : dict
        Column name (must be in 'metadata') to palette name mapping
    color_cols : list
        Column names in 'metadata' to color by

    Returns
    -------
    graph_dict : dict of holoviews.Graph
        (alphabet, ksize, log2sketchsize) : holoviews.Graph mapping for all
        similarity matrices found. To be used by 'draw_holoviews_graphs'

    """
    # Strip the final slash because it makes s3 stuff weird
    figure_folder = figure_folder.rstrip("/")

    iterable = itertools.product(molecules, ksizes, log2sketchsizes)
    graph_dict = {}

    categories = metadata[color_cols]

    for molecule, ksize, log2sketchsize in iterable:
        template_kwargs = dict(
            molecule=molecule, ksize=ksize, log2sketchsize=log2sketchsize
        )
        sketch_id = sketch_id_template.format(**template_kwargs)
        if verbose:
            print(sketch_id.replace("-", ": ").replace("_", ", "))
        csv = csv_template.format(**template_kwargs)
        try:
            similarities = pd.read_csv(csv)
        except FileNotFoundError:
            warnings.warn(f"file {csv} not found")
            # File doesn't exist yet
            continue
        similarities.index = similarities.columns
        if verbose:
            print(f"\tsimilarities.shape: {similarities.shape}")

        title = (
            f"alphabet: {molecule}, ksize: {ksize}, "
            f"log2sketchsize: {log2sketchsize}"
        )

        if plaidplot:
            try:
                g = sourmash_utils.plaidplot(
                    similarities,
                    metric="cosine",
                    row_categories=categories,
                    col_categories=categories,
                    row_palette=palettes,
                    col_palette=palettes,
                )
                g.fig.suptitle(title)
                png = f"{figure_folder}/{sketch_id}_plaidplot.png"
                savefig(g, png, dpi=150)
            except FloatingPointError:
                warnings.warn("\tCouldn't compute linkage -- no plaidplot " "generated")

        graph, pos = build_graph_and_plot(
            similarities,
            metadata,
            n_neighbors,
            color_cols,
            palettes,
            figure_folder,
            sketch_id,
            title,
        )

        # hv.extension('matplotlib')

        graph_hv = hv.Graph.from_networkx(graph, pos)

        graph_hv = graph_hv.opts(
            node_size=10,
            edge_line_width=1,
            cmap="Set2",
            node_color=dim(groupby),
            node_line_color="gray",
        )
        bundled = bundle_graph(graph_hv)
        # hv.save(bundled, '.pdf', backend='matplotlib')
        graph_dict[(molecule, ksize, log2sketchsize)] = bundled

        if make_within_groupby_graphs:
            # make within-group (e.g. within-species) graphs
            for species, df in metadata.groupby(groupby):
                data = similarities.loc[df.index, df.index]
                figure_prefix = f"{sketch_id}_{species}"
                graph_title = f"{title} ({species})"
                build_graph_and_plot(
                    data,
                    df,
                    n_neighbors,
                    color_cols,
                    palettes,
                    figure_folder,
                    figure_prefix,
                    graph_title,
                )

    return graph_dict


def draw_holoviews_graphs(graph_dict):
    # use first key to determine default settings
    first_key = list(graph_dict.keys())[0]
    molecule, ksize, log2sketchsize = first_key

    hv.extension("bokeh")
    defaults = dict(width=400, height=400, padding=0.1)
    hv.opts.defaults(
        opts.EdgePaths(**defaults), opts.Graph(**defaults), opts.Nodes(**defaults)
    )

    kdims = [
        hv.Dimension(("alphabet", "alphabet"), default=molecule),
        hv.Dimension(("ksize", "k-mer size"), default=ksize),
        hv.Dimension(("log2_num_hashes", "$log_2$ num hashes"), default=log2sketchsize),
    ]

    kwargs = dict(width=800, height=800, xaxis=None, yaxis=None)
    opts.defaults(opts.Nodes(**kwargs), opts.Graph(**kwargs))

    kwargs = dict(
        node_size=10,
        edge_line_width=1,
        cmap="Set2",
        node_color=dim("species"),
        node_line_color="gray",
        width=600,
        height=600,
        xaxis=None,
        yaxis=None,
    )

    holomap = hv.HoloMap(graph_dict, kdims=kdims)
    holomap.opts(opts.Graph(**kwargs))
    return holomap
