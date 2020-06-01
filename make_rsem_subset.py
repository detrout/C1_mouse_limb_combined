#!/usr/bin/env python3
from argparse import ArgumentParser
import logging
import os
import pandas

from woldrnaseq import models
from woldrnaseq.makersemcsv import GeneRsemLoader

from generate_combined_transcript_C1 import (
    ASOF_RUN17_library_files,
    ASOF_RUN17_experiment_files
)


logger = logging.getLogger('make_rsem_subset')


def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('filename',
                        help='name of cell set spreadsheet to load')
    parser.add_argument('--gtf-cache', help='Specify name of GTF cache file')
    parser.add_argument('--add-names', action='store_true', default=False,
                        help='Add names to ouptut quantification file')
    args = parser.parse_args()

    if args.add_names:
        if args.gtf_cache is None:
            parser.error('GTF-cache is needed to add names to the quantification file')
        else:
            logger.info('Loading GTF Cache %s', args.gtf_cache)
            annotation = models.load_gtf_cache(args.gtf_cache)
    else:
        annotation = None

    libraries = load_asof_run17_libraries()
    subset = load_cells_set(args.filename)
    experiment = pandas.Series({
        'replicates': list(subset.index)})
    experiment.name = 'subset'

    print(experiment.head())
    print(experiment.replicates[0], type(experiment.replicates[0]))
    loader_tpm = GeneRsemLoader('TPM', annotation=annotation)
    loader_length = GeneRsemLoader('length', annotation=annotation)
    loader_effective_length = GeneRsemLoader('effective_length', annotation=annotation)
    matrix_tpm = loader_tpm.load(experiment, libraries)
    matrix_length = loader_length.load(experiment, libraries)
    matrix_effective_length = loader_effective_length.load(experiment, libraries)

    if args.add_names:
        insert_location = 1
    else:
        insert_location = 0
    length = matrix_length[matrix_length.columns[-1]]
    matrix_tpm.insert(insert_location, 'gene_length', length)
    print(matrix_tpm.head())
    matrix_tpm.to_csv('transcript_length_mouse_pacbio_pools.TPM.csv')
    matrix_effective_length.to_csv('transcript_length_mouse_pacbio_pools.effective_length.csv')


def split_files_text(text):
    for row in text.split('\n'):
        yield os.path.expanduser(row.strip())


def sanitize_library_name(name):
    return name.replace('_mm10', '').replace('_clean', '')


def load_asof_run17_libraries():
    library_files = list(split_files_text(ASOF_RUN17_library_files))
    libraries = models.load_library_tables(library_files)
    name = libraries.index.name
    libraries.index = [sanitize_library_name(x) for x in libraries.index]
    libraries.index.name = name

    return libraries


def load_asof_run17_experiments():
    experiment_files = list(split_files_text(ASOF_RUN17_experiment_files))
    experiments = models.load_experiments(experiment_files)

    for experiment_name, row in experiments.iterrows():
        row.replicates = [sanitize_library_name(x) for x in row.replicates]

    return experiments


def load_cells_set(filename, sheet=0):
    if filename.endswith('.xlsx'):
        df = pandas.read_excel(filename, sheet=sheet)
    elif filename.endswith('.tsv'):
        df = pandas.read_csv(filename, sep='\t')
    else:
        raise ValueError('Unrecognized filename extension {}'.format(filename))

    df.columns = ['library_id', 'cluster_assignment']
    df.set_index('library_id', inplace=True)
    filtered = df.dropna(axis=0, how='any')
    print('Filtered to cluster members', filtered.shape)
    return filtered


if __name__ == "__main__":
    main()

