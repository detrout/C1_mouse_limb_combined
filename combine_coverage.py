"""Combine coverage files for all of our C1 mouse limb data
"""
from glob import glob
from pprint import pprint
import pandas
import os
import sys

from woldrnaseq.models import load_coverage

from to_include import (
    generate_to_include,
    generate_to_include_as_of_run17
)

roots_433 = """~sau/flowcells/C1_e10.5_mouse_limb_run2_June20_2016/
~sau/flowcells/H5LV3BCXY/
~diane/proj/C1_e10.5_mouse_limb_run3_Dec5_2016_2/
~sau/flowcells/C1_mouse_limb_combined_Mar_2017/
"""

roots_asof_run17 = roots_433 + """~sau/flowcells/HFNLNBCX2/
~sau/flowcells/H7CNTBCX2/
~sau/flowcells/HFNLTBCX2/
~sau/flowcells/HF7NTBCX2/
~sau/flowcells/HFNYNBCX2/
~sau/flowcells/HJFCLBCXY/
"""

def main():
    make_asof_run17_coverage()


def make_433_coverage():
    # to_includes includes gene_id header, remove it. with [1:]
    coverage = load_all_coverage(roots_433, generate_to_include()[1:])
    #coverage_df.to_csv('C1_mouse_limb_combined_coverage.tsv', sep='\t')


def make_asof_run17_coverage():
    # to_includes includes gene_id header, remove it. with [1:]
    coverage = load_all_coverage(roots_asof_run17, generate_to_include_as_of_run17()[1:])
    coverage.to_csv("C1_mouse_limb_coverage_asof_run17.tsv", sep='\t')


def load_all_coverage(roots, to_include):
    sheets = None
    coverage = []
    found = set()
    for path in roots.split():
        for lib in to_include:
            coverage_pattern = os.path.join(os.path.expanduser(path), lib, '*.coverage')
            coverage_files = glob(coverage_pattern)
            if len(coverage_files) == 1:
                coverage.append(load_coverage(coverage_files[0], lib))
                #print(lib, coverage_files[0])
                found.add(lib)
            elif len(coverage_files) > 1:
                print('too many for {}'.format(lib))
            #print(path.split('/')[-2], lib, len(coverage))

    print('not found', set(to_include).difference(found))
    coverage_df = pandas.concat(coverage, axis=1)
    print(coverage_df.shape)
    return coverage_df

if __name__ == "__main__":
    main()
