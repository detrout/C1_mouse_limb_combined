"""

"""
import seaborn
from matplotlib import pyplot as plt
import requests
from io import BytesIO
import os
import sys


PANDASODF = os.path.expanduser('~diane/src/pandasodf')
if PANDASODF not in sys.path:
    sys.path.append(PANDASODF)
from pandasodf import ODFReader


class C1BoxPlotter(seaborn.categorical._BoxPlotter):
    """Custom seaborn plot for C1 box plot figures"""
    def restyle_boxplot(self, artist_dict, color, props):
        super(C1BoxPlotter, self).restyle_boxplot(artist_dict, color, props)
        color = '#C0C0C0'
        if artist_dict['boxes'][0].get_facecolor() == (0.0, 0.0, 0.0, 1.0):
            artist_dict['medians'][0].set_color(color)


def boxplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
            orient=None, color=None, palette=None, saturation=.75,
            width=.8, dodge=True, fliersize=5, linewidth=None,
            whis=1.5, notch=False, ax=None, **kwargs):
    """Custom seaborn plot for C1 box plot figures

    They want one of the clusters to be black, so we need to
    turn the median line white.
    """

    plotter = C1BoxPlotter(x, y, hue, data, order, hue_order,
                           orient, color, palette, saturation,
                           width, dodge, fliersize, linewidth)

    if ax is None:
        ax = plt.gca()
    kwargs.update(dict(whis=whis, notch=notch))

    plotter.plot(ax, kwargs)
    return ax


def read_remote_sheet(url, sheet_name):
    resp = requests.get(url)
    content = BytesIO(resp.content)
    book = ODFReader(content)
    sheet = book.parse(
        sheetname=sheet_name,
        dtype={
            '10x_class': int,
            'order': int,
            'Red': int,
            'Green': int,
            'Blue': int,
        })
    return sheet


def get_cluster_maps(sheet, class_column, label_column):
    filtered_order = sheet.sort_values('order')[[class_column, label_column]].dropna()
    return {
        'label': dict(sheet[[class_column, label_column]].dropna().values),
        'color': dict(sheet[[label_column, 'color']].dropna().values),
        'order': filtered_order[label_column].values,
    }
