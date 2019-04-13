import itertools
import warnings

# 3rd-party modules
import holoviews as hv
from holoviews import opts, dim
from holoviews.operation.datashader import datashade, bundle_graph
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import seaborn as sns
%matplotlib inline

# My handwritten modules
from aguamenti.s3_utils import write_s3
from . import knn
from . import sourmash_utils
from . import sourmash_compare_utils


hv.extension('bokeh')

defaults = dict(width=400, height=400, padding=0.1)
hv.opts.defaults(
    opts.EdgePaths(**defaults), opts.Graph(**defaults), opts.Nodes(**defaults))

import common

figure_folder = common.get_figure_folder()
! mkdir -p $figure_folder
figure_folder

# don't warn me about too many figures open
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})



KSIZES = 9, 12, 15, 21
LOG2SKETCHSIZES = 10, 12, 14, 16
MOLECULES = 'dna', 'protein'

color_cols = ['species', 'cell_label', ]
categories = metadata_source_name[color_cols]

palettes = dict(species='Set2', cell_label='tab20')

s3_folder = 's3://olgabot-maca/nf-kmer-similarity/human_mouse_zebrafish/'
SKETCH_ID_TEMPLATE = 'molecule-{molecule}_ksize-{ksize}_log2sketchsize-{log2sketchsize}'
template = s3_folder + f
'similarities_{sketch_id_template}.csv'

# bundled_graphs = []
# graph_dict = {}

n_neighbors = 5


def build_graph_and_plot(data, metadata, n_neighbors, color_cols, palettes,
                         figure_folder, figure_prefix, title):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        graph = knn.nearest_neighbor_graph(data, metadata,
                                           n_neighbors=n_neighbors,
                                           color_cols=color_cols,
                                           palettes=palettes)

    pos = nx.spring_layout(graph, seed=0)

    for label in color_cols:
        fig, ax = plt.subplots()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            knn.draw_graph(graph, edge_color='black', label_col=label, pos=pos)
        ax.set_title(title)
        png = f
        '{figure_folder}/{figure_prefix}_graph_nneighbors-{n_neighbors}_colorby-{label}.png'
        fig.savefig(png, dpi=150)


def read_similarity_matrices(csv_template, ksizes=KSIZES,
                             log2sketchsizes=LOG2SKETCHSIZES,
                             molecules=MOLECULES,
                             sketch_id_template=SKETCH_ID_TEMPLATE):
    iterable = itertools.product(molecules, ksizes, log2sketchsizes)

    for ksize, log2sketchsize in iterable:
        template_kwargs = dict(ksize=ksize, log2sketchsize=log2sketchsize)
        sketch_id = sketch_id_template.format(**template_kwargs)
        print(sketch_id.replace('-', ": ").replace("_", ", "))
        csv = csv_template.format(**template_kwargs)
        try:
            similarities = pd.read_csv(csv)
        except FileNotFoundError:
            print(f"file {csv} not found")
            # File doesn't exist yet
            continue
        similarities.index = similarities.columns
        print(similarities.shape)

        title = f
        "ksize: {ksize}, log2sketchsize: {log2sketchsize}"

        palette = dict(species='tab10', cell_ontology_class='husl')

        #     g = sourmash_utils.plaidplot(similarities,
        #                                  metric='cosine', row_categories=categories,
        #                              col_categories=categories, row_palette=palettes, col_palette=palettes)
        #     g.fig.suptitle(title)
        #     png = f'{figure_folder}/{sketch_id}_plaidplot.png'
        #     g.savefig(png, dpi=150)

        build_graph_and_plot(similarities, metadata_source_name,
                             n_neighbors, color_cols, palettes, figure_folder,
                             sketch_id, title)

        #     graph = knn.nearest_neighbor_graph(similarities, metadata_source_name, n_neighbors=n_neighbors,
        #                                       color_cols=color_cols, palettes=palettes)

        #     pos = nx.spring_layout(graph, seed=0)

        #     for label in color_cols:
        #         fig, ax = plt.subplots()
        #         knn.draw_graph(graph, edge_color='black', label_col=label, pos=pos)
        #         ax.set_title(title)
        #         png = f'{figure_folder}/{sketch_id}_graph_nneighbors-{n_neighbors}_colorby-{label}.png'
        #         fig.savefig(png, dpi=150)

        # make within-species graphs
        for species, df in metadata_source_name.groupby('species'):
            data = similarities.loc[df.index, df.index]
            build_graph_and_plot(data, df,
                                 n_neighbors, color_cols, palettes, figure_folder,
                                 f
            "{sketch_id}_{species}", title = f
            "{title} ({species})")

            #     graph_hv = hv.Graph.from_networkx(graph, pos)

            #     graph_hv = graph_hv.opts(node_size=10, edge_line_width=1, cmap='Set2', node_color=dim("species"),
            #                           node_line_color='gray')
            #     bundled = bundle_graph(graph_hv)
            #     hv.save(bundled, 'asdf.pdf', backend='matplotlib')
            #     bundled_graphs.append(bundled)
            #     graph_dict[(ksize, log2sketchsize)] = bundled

            # graph_dict = dict(zip(iterable, bundled_graphs))



