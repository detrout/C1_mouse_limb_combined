"""Combine coverage files for all of our C1 mouse limb data
"""
from glob import glob
from pprint import pprint
import pandas
import os
import sys

from woldrnaseq.models import load_coverage


libraries = [
    '15019',
    '15020',
    '15084',
    '15085',
    '16111',
    '16112',
    '16134',
    '16135',
    '16930',
    '16931',
    '17298',
    '17299',
]
roots = """/woldlab/castor/home/diane/proj/C1_mouse_limb_combined/bulk/
"""

sheets = None
coverage = []
for path in roots.split():
    for lib in libraries:
        coverage_pattern = os.path.join(path, lib, '*.coverage')
        coverage_files = glob(coverage_pattern)
        if len(coverage_files) == 1:
            coverage.append(load_coverage(coverage_files[0], lib))
        elif len(coverage_files) > 1:
            print('too many for {}'.format(lib))
        print(path.split('/')[-2], lib, len(coverage))

coverage_df = pandas.concat(coverage, axis=1)
print(coverage_df.shape)
#pprint(sorted(coverage_df.columns))
coverage_df.to_csv('bulk_mouse_limb_combined_coverage.tsv', sep='\t')
