#!/usr/bin/python3
"""Link data files from scattered analysis
"""
from argparse import ArgumentParser
import os

from woldrnaseq.models import (
    load_experiments,
    load_library_tables,
    genome_name_from_library,
)
from woldrnaseq.make_tracks import (
    make_bam_track_name,
)
from generate_combined_transcript_C1 import ASOF_RUN17_experiment_files, ASOF_RUN17_library_files
from find_bigwigs import (
    read_peng_20180710_cluster_memberships,
    sanitize_library_suffix,
)

def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('-o', '--output-dir')
    args = parser.parse_args(cmdline)

    experiment_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_experiment_files.split() ]
    library_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_library_files.split() ]

    experiments = load_experiments(experiment_files)
    libraries = load_library_tables(library_files)

    #link_rsem(libraries, args.output_dir)
    link_genome_bams(libraries, args.output_dir)


def link_rsem(libraries, output_dir):
    for library_id, library in libraries.iterrows():
        clean_library_id = sanitize_library_suffix(library_id)
        target_dir = os.path.join(output_dir, clean_library_id)
        if not os.path.exists(target_dir):
            os.mkdir(target_dir)

        source_dir = library.analysis_dir
        cur_dir = os.getcwd()
        os.chdir(target_dir)
        for extension in ['_anno_rsem.genes.results', '_anno_rsem.isoforms.results']:
            suffix = '-' + genome_name_from_library(library) + extension
            source_name = library_id + suffix
            target_name = clean_library_id + suffix
            source_pathname = os.path.join(source_dir, source_name)
            if os.path.exists(source_pathname) and not os.path.exists(target_name):
                print(source_pathname, '->', target_name)
                os.symlink(source_pathname, target_name)
        os.chdir(cur_dir)

def link_genome_bams(libraries, output_dir):
    for library_id, library in libraries.iterrows():
        clean_library_id = sanitize_library_suffix(library_id)
        target_dir = os.path.join(output_dir, clean_library_id)
        if not os.path.exists(target_dir):
            os.mkdir(target_dir)

        name = make_bam_track_name(library, library.analysis_dir)
        source_pathname = os.path.join(library.analysis_dir, name)
        target_name = clean_library_id + genome_name_from_library(library) + '_genome.bam'
        cur_dir = os.getcwd()
        os.chdir(target_dir)
        if os.path.exists(source_pathname) and not os.path.exists(target_name):
            print(source_pathname, '->', target_name)
            os.symlink(source_pathname, target_name)
        os.chdir(cur_dir)


if __name__ == '__main__':
    main()
