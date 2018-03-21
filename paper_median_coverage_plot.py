"""Generate plots median plots brian asked for for paper

One showing all 433 cells, and one for bulk data
"""
import pandas

from woldrnaseq.plot_coverage import make_median_normalized_summary

def main(cmdline=None):
    make_generic_plot('C1_mouse_limb_combined_coverage.tsv', 'C1_mouse_limb')
    make_generic_plot('bulk_mouse_limb_combined_coverage.tsv', 'bulk_mouse_limb')


def make_generic_plot(filename, experiment_name):
    table = pandas.read_csv(filename,
                             sep='\t',
                             index_col=0,
    )
    experiment = { experiment_name: list(table.keys()) }

    print(table.shape)
    print(table.describe())
    plot = make_median_normalized_summary(experiment, table)

if __name__ == '__main__':
    main()
