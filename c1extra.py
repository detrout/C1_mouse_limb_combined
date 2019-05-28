"""Custom seaborn plot for C1 box plot figures

They want one of the clusters to be black, so we need to
turn the median line white.
"""
import seaborn
from matplotlib import pyplot as plt


class C1BoxPlotter(seaborn.categorical._BoxPlotter):
    def restyle_boxplot(self, artist_dict, color, props):
        super(C1BoxPlotter, self).restyle_boxplot(artist_dict, color, props)
        color = '#C0C0C0'  #'#717E8D'
        if artist_dict['boxes'][0].get_facecolor() == (0.0, 0.0, 0.0, 1.0):
            #artist_dict['boxes'][0].set_edgecolor(color)
            artist_dict['medians'][0].set_color(color)


def boxplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
            orient=None, color=None, palette=None, saturation=.75,
            width=.8, dodge=True, fliersize=5, linewidth=None,
            whis=1.5, notch=False, ax=None, **kwargs):

    plotter = C1BoxPlotter(x, y, hue, data, order, hue_order,
                           orient, color, palette, saturation,
                           width, dodge, fliersize, linewidth)

    if ax is None:
        ax = plt.gca()
    kwargs.update(dict(whis=whis, notch=notch))

    plotter.plot(ax, kwargs)
    return ax
