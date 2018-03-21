import pandas
import os
import sys

from to_include import generate_to_include

to_include = generate_to_include()

files = """/woldlab/castor/home/sau/flowcells/C1_e10.5_mouse_limb_run2_June20_2016/C1_e10.5_mouse_limb_run2_June20_2016_FPKM.csv
/woldlab/castor/home/sau/flowcells/H5LV3BCXY/C1_e10.5_mouse_limb_run1_June6_2016_FPKM.csv
/woldlab/castor/home/diane/proj/C1_e10.5_mouse_limb_run3_Dec5_2016_2/C1_e10.5_mouse_limb_run3_Dec5_2016_gene_FPKM.csv
/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/C1_mouse_e13.5_limb_mesenchyme_mm10_run4_gene_FPKM.csv
/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/C1_mouse_e11.0_limb_mesenchyme_mm10_run5_gene_FPKM.csv
/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/C1_mouse_e11.5_limb_mesenchyme_mm10_run6_gene_FPKM.csv
/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/C1_mouse_e12.5_limb_mesenchyme_mm10_run7_gene_FPKM.csv
/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017/C1_mouse_e13.5_limb_mesenchyme_mm10_run8_gene_FPKM.csv"""


sheets = None
for f in files.split():
    if sheets is None:
        sheets = pandas.read_csv(f)
    else:
        sheets = pandas.merge(sheets, pandas.read_csv(f), on='gene_id')

print(sheets.columns)
print(sheets.shape)
print(sheets[to_include].shape)
filtered = sheets[to_include]
newsheet = filtered.set_index('gene_id')

newsheet.to_csv('C1_mouse_combined_peng_filter.csv')
