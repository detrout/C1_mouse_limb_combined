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

from combine_coverage import (
    roots_433,
    roots_asof_run17
)

def main():
    detected_files = find_detected_files()
    if detected_files is not None:
        print(detected_files.shape)
        #print(detected_files)
        detected_files.to_csv("genes-detected-C1_mouse_limb_combined_asof_run17.tsv", sep='\t')

def find_detected_files():
    ignore = [
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e10.5_limb_mm10_clean_run1_june_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e10.5_limb_mm10_clean_run1_seed_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_e10.5_mouse_limb_run2_June20_2016/genes-detected_C1_e10.5_mouse_limb_run2_June20_2016_FPKM.csv',
#        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e10.5_limb_mm10_clean_run2_gene_FPKM.csv',
        '/home/diane/proj/C1_e10.5_mouse_limb_run3_Dec5_2016_2/genes-detected_C1_e10.5_mouse_limb_run3_Dec5_2016_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e13.5_limb_mesenchyme_run4_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e13.5_limb_mesenchyme_mm10_clean_run4_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e13.5_limb_mesenchyme_run4_gene_FPKM.csv'
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e13.5_limb_mesenchyme_run4_gene_FPKM.csv'
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e13.5_limb_mesenchyme_mm10_clean_run4_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e11.0_limb_mesenchyme_run5_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e11.0_limb_mesenchyme_mm10_clean_run5_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e11.5_limb_mesenchyme_mm10_clean_run6_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e11.5_limb_mesenchyme_run6_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e12.5_limb_mesenchyme_run7_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e12.5_limb_mesenchyme_mm10_clean_run7_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e13.5_limb_mesenchyme_mm10_clean_run8_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/genes-detected_C1_mouse_e13.5_limb_mesenchyme_run8_gene_FPKM.csv',
        '/woldlab/castor/home/sau/flowcells/HJFCLBCXY/genes-detected_C1_mouse_e12.5_limb_run9_March8_2017_poolC_HJFCLBCXY_gene_FPKM.csv',
    ]
    seen = {}
    tables = []
    for d in roots_asof_run17.split():
        paths = os.path.expanduser(os.path.join(d, 'genes-detected*'))
        for detected in glob(paths):
            if detected in ignore:
                #print('skip')
                continue

            print(detected)
            try:
                table = pandas.read_csv(detected, index_col=0)
                library_ids = []
                for library_id in table.index:
                    library_id = library_id.replace('_mm10', '').replace('_clean', '').replace('_june', '').replace('_seed', '')
                    library_ids.append(library_id)
                    seen.setdefault(library_id, []).append(detected)

                table.index = library_ids
                tables.append(table)
            except pandas.errors.ParserError:
                print('unable to read', detected)

    dups = 0
    for key in sorted(seen):
        if len(seen[key]) > 1:
            print(key, seen[key])
            dups += 1

    if dups == 0:
        return pandas.concat(tables)


if __name__ == '__main__':
    main()
