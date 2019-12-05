from matplotlib.pyplot import plt
import seaborn as sns
from umap import UMAP


def plot_umap(transformed, *args, **kwargs):
    x = transformed[:, 0]
    y = transformed[:, 1]

    fig, ax = plt.subplots()
    sns.scatterplot(x, y, *args, **kwargs)
    sns.despine(left=True, bottom=True)
    ax.set(xlabel='UMAP 1', ylabel='UMAP 2', xticks=[], yticks=[])
    return ax


def do_umap(similarities, random_state=0, y=None):
    umapper = UMAP(metric='precomputed', random_state=random_state,
                   *args, **kwargs)
    transformed = umapper.fit_transform(1-similarities, y=y)
    return transformed


def umap_and_plot(similarities, random_state=0, hue=None, *args, **kwargs):
    transformed = do_umap(similarities, random_state=random_state)
    plot_umap(transformed, hue=hue, *args, **kwargs)
