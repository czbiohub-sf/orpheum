import os


import jupyter_utils



FIGURE_FOLDER = os.path.join('..', 'figures')

def get_notebook_basename():
    notebook_path = jupyter_utils.get_notebook_name()
    notebook_basename = os.path.basename(notebook_path).split('.ipynb')[0]
    return notebook_basename


def get_figure_folder():
    notebook_basename = get_notebook_basename()
    figure_folder = os.path.join(FIGURE_FOLDER, notebook_basename)
    return figure_folder