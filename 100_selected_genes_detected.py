"""Generate plots median plots brian asked for for paper

This 
"""
import pandas

from woldrnaseq.plot_coverage import make_median_normalized_summary

def main(cmdline=None):
    list_of_100 = pandas.read_excel('list_of_100_single_cell_libraries_for_genes_called.xlsx')
    list_of_100 = [ x[1:] for x in list_of_100['library #']]

    table = pandas.read_csv('C1_mouse_combined_peng_filter.tsv',
                             sep='\t',
                             index_col=0,
    )
    print(table.columns)
    print(table.index)
    table = table[list_of_100]

    table.to_csv('list_of_100_single_cell_libraries.csv', )

if __name__ == '__main__':
    main()
