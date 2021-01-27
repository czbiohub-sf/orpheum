import os

import pandas as pd

from . import jupyter_utils


FIGURE_FOLDER = os.path.join("..", "figures")


def get_notebook_basename():
    notebook_path = jupyter_utils.get_notebook_name()
    notebook_basename = os.path.basename(notebook_path).split(".ipynb")[0]
    return notebook_basename


def get_figure_folder():
    notebook_basename = get_notebook_basename()
    figure_folder = os.path.join(FIGURE_FOLDER, notebook_basename)
    return figure_folder


def load_tabula_muris_annotations():
    """Fetch latest FACS annotations of Tabula Muris data from GitHub"""
    annotations = pd.read_csv(
        "https://github.com/czbiohub/tabula-muris/raw/master/00_data_ingest/"
        "18_global_annotation_csv/annotations_facs.csv",
        index_col="cell",
    )

    # Sanitize the input - no dots
    annotations.index = annotations.index.str.replace(".", "-")
    annotations.columns = annotations.columns.str.replace(".", "_")
    annotations["sample_id"] = annotations.index
    annotations = annotations.fillna("NA")
    return annotations
