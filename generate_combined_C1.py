import pandas
import os
import sys

from to_include import generate_to_include, generate_to_include_as_of_run17


paper_433_FPKM_files = """~sau/flowcells/C1_e10.5_mouse_limb_run2_June20_2016/C1_e10.5_mouse_limb_run2_June20_2016_FPKM.csv
    ~sau/flowcells/H5LV3BCXY/C1_e10.5_mouse_limb_run1_June6_2016_FPKM.csv
    ~diane/proj/C1_e10.5_mouse_limb_run3_Dec5_2016_2/C1_e10.5_mouse_limb_run3_Dec5_2016_gene_FPKM.csv
    ~sau/flowcells/C1_mouse_limb_combined_Mar_2017/C1_mouse_e13.5_limb_mesenchyme_mm10_run4_gene_FPKM.csv
    ~sau/flowcells/C1_mouse_limb_combined_Mar_2017/C1_mouse_e11.0_limb_mesenchyme_mm10_run5_gene_FPKM.csv
    ~sau/flowcells/C1_mouse_limb_combined_Mar_2017/C1_mouse_e11.5_limb_mesenchyme_mm10_run6_gene_FPKM.csv
    ~sau/flowcells/C1_mouse_limb_combined_Mar_2017/C1_mouse_e12.5_limb_mesenchyme_mm10_run7_gene_FPKM.csv
    ~sau/flowcells/C1_mouse_limb_combined_Mar_2017/C1_mouse_e13.5_limb_mesenchyme_mm10_run8_gene_FPKM.csv"""

ASOF_RUN17_FPKM_FILES = paper_433_FPKM_files + """
    ~sau/flowcells/HFNLNBCX2/C1_mouse_e12.5_forelimb_run9_March8_2017_rearray_March152018_gene_FPKM.csv
    ~sau/flowcells/H7CNTBCX2/C1_mouse_e11.5_forelimb_run10_December11_2017_gene_FPKM_names.csv
    ~sau/flowcells/HFNLTBCX2/C1_mouse_e12.0_forelimb_run11_December12_2017_gene_FPKM_names.csv
    ~sau/flowcells/HF7NTBCX2/C1_mouse_e12.0_forelimb_run11_December12_2017_gene_FPKM_names.csv
    ~sau/flowcells/HF7NTBCX2/C1_mouse_e13.0_forelimb_run12_December13_2017_gene_FPKM_names.csv
    ~sau/flowcells/HFNLNBCX2/C1_mouse_e14.0_forelimb_run13_December14_2017_gene_FPKM_names.csv
    ~sau/flowcells/HFNYNBCX2/C1_mouse_e15.5_forelimb_run14_December15_2017_gene_FPKM_names.csv
    ~sau/flowcells/H7CNTBCX2/C1_mouse_e10.5_forelimb_run15_January13AM_2018_gene_FPKM_names.csv
    ~sau/flowcells/H7CNTBCX2/C1_mouse_e11.0_forelimb_run16_January13PM_2018_gene_FPKM_names.csv
    ~sau/flowcells/H7CNTBCX2/C1_mouse_e14.5_forelimb_run17_January16_2018_gene_FPKM_names.csv
"""

def main():
    #make_paper_combined()
    make_asof_run17_combined()


def make_paper_combined():
    files = paper_433_FPKM_files
    to_include = generate_to_include()

    combined = make_combined(files, to_include)
    combined.to_csv('C1_mouse_combined_peng_filter.tsv', sep='\t')


def make_asof_run17_combined():
    files = ASOF_RUN17_FPKM_FILES
    to_include = generate_to_include_as_of_run17()

    combined = make_combined(files, to_include)
    combined.to_csv('C1_mouse_combined_asof_run17.tsv', sep='\t')


def make_combined(quant_files, to_include):
    sheets = None
    gene_name = None
    for f in quant_files.split():
        f = os.path.expanduser(f.strip())
        if sheets is None:
            sheets = pandas.read_csv(f)
        else:
            sheet = pandas.read_csv(f)
            columns = [ x.replace('_mm10', '').replace('_clean', '') for x in sheet.columns]
            sheet.columns = columns
            if 'gene_name' in sheet.columns:
                gene_name = sheet.gene_name
                sheet = sheet.drop('gene_name', axis=1)
            sheets = pandas.merge(sheets, sheet, on='gene_id')

    print(sheets.columns)
    print(sheets.shape)
    print('to_include', to_include)
    print('columns', list(sheets.columns))
    print(sheets[to_include].shape)
    filtered = sheets[to_include]
    filtered = sheets
    newsheet = filtered.set_index('gene_id')
    print(newsheet.shape)
    return newsheet

if __name__ == '__main__':
    main()
