import pandas
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot


def main():
    coverage = pandas.read_csv('C1_mouse_limb_coverage_asof_run17.tsv', sep='\t')

    slopes = {}
    for i, column in enumerate(coverage):
        line = coverage[column]
        slope = compute_coverage_slope(line)
        slopes[column] = slope
        print(i, len(coverage.columns), column, slope)

    slopes = pandas.Series(slopes)
    slopes.to_csv('asof_run17_coverage_slopes.tsv', sep='\t')
    #print(slopes.describe())
    #print(numpy.histogram(slopes.values))
    f = pyplot.figure()
    slopes.hist(ax=f.gca())
    f.savefig('asof_run17_slopes_hist.png')
    
def compute_coverage_slope(line):
    plateau = line[20:81]
    m, b = numpy.polyfit(numpy.arange(20, 81), plateau, 1)
    return m
    
if __name__ == '__main__':
    main()
