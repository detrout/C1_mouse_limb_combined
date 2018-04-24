import pandas
import os
import sys

from woldrnaseq import models
from woldrnaseq.makersemcsv import IsoformRsemLoader

from to_include import generate_to_include, generate_to_include_asof_run17

paper_433_library_files = """~sau/flowcells/H5LV3BCXY/library.tsv
    ~sau/flowcells/C1_e10.5_mouse_limb_run2_June20_2016/library.tsv
    ~diane/proj/C1_e10.5_mouse_limb_run3_Dec5_2016_2/library.tsv
    ~sau/flowcells/C1_mouse_limb_combined_Mar_2017/library_mm10.tsv"""

ASOF_RUN17_library_files = paper_433_library_files + """
    ~sau/flowcells/HFNLNBCX2/library.tsv
    ~sau/flowcells/H7CNTBCX2/library.tsv
    ~sau/flowcells/HFNLTBCX2/library.tsv
    ~sau/flowcells/HF7NTBCX2/library.tsv
    ~sau/flowcells/HFNYNBCX2/library.tsv"""

paper_433_experiment_files = """~sau/flowcells/H5LV3BCXY/experiments.tsv
    ~sau/flowcells/C1_e10.5_mouse_limb_run2_June20_2016/all-in-one-experiment.tsv
    ~diane/proj/C1_e10.5_mouse_limb_run3_Dec5_2016_2/experiments.tsv
    ~sau/flowcells/C1_mouse_limb_combined_Mar_2017/experiments_mm10.tsv"""

ASOF_RUN17_experiment_files = paper_433_experiment_files + """
    ~sau/flowcells/HFNLNBCX2/experiments.tsv
    ~sau/flowcells/H7CNTBCX2/experiments.tsv
    ~sau/flowcells/HFNLTBCX2/experiments.tsv
    ~sau/flowcells/HF7NTBCX2/experiments.tsv
    ~sau/flowcells/HFNYNBCX2/experiments.tsv"""


def main():
    print('loading gtf')
    gtf = models.load_gtf_cache('/woldlab/castor/home/sau/genomes/mm10-M4-male/mm10-M4-male.h5')
    print('loading isoforms')
    filtered = load_filtered_transcripts()
    # toi = transcripts of interest
    toi = gtf[[gigios_genes_filter(x) for x in gtf['gene_id']]]['transcript_id'].values
    transcripts = filtered.loc[toi]
    print('saving', transcripts.shape)
    
    transcripts.to_csv('gigios_transcripts_fpkms.tsv', sep='\t')

def load_filtered_transcripts():
    sep = '\t'

    cache_file = os.path.expanduser('~sau/genomes/mm10-M4-male/mm10-M4-male.h5')
    #annotation = models.load_gtf_cache(cache_file)
    annotation = None

    loader = IsoformRsemLoader('FPKM', annotation)
    index_name = 'transcript_id'
    # loader = GeneRsemLoader(args.quantification, annotation)
    #index_name = 'gene_id'
    
    to_include = generate_to_include_asof_run17()[1:]

    experiment_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_experiment_files.split() ]
    library_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_library_files.split() ]
    
    quantifications = []
    for e, l in zip(experiment_files, library_files):
        print('loading', e)
        experiments = models.load_experiments([e], sep=sep)
        libraries = models.load_library_tables([l], sep=sep)
        for i, experiment in experiments.iterrows():
            print(experiment)
            quantification = loader.load(experiment, libraries)
            quantification.columns = list(filter_columns(quantification.columns))
            quantifications.append(quantification)

    sheets = pandas.concat(quantifications, axis=1)

    print('all', sheets.shape)
    # sheets.to_csv('C1_mouse_combined_transcript_asof_run17_unfiltred.tsv', sep='\t')
    # was crashing because of _mm10 suffix
    filtered = sheets[to_include]
    print('filtered', filtered.shape)
    return filtered
    #filtered.to_csv('C1_mouse_combined_transcript_asof_run17.tsv', sep='\t')
    

def filter_columns(columns):
    for c in columns:
        c = c.replace('_mm10', '').replace('_clean', '')
        yield c

def gigios_genes_filter(x):
    for gene_id in ['ENSMUSG00000020167','ENSMUSG00000032228','ENSMUSG00000063659','ENSMUSG00000030189']:
        if x.startswith(gene_id):
            return True
    return False

if __name__ == '__main__':
    main()
