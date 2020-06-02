#!/usr/bin/python3
from argparse import ArgumentParser
import os
import logging
import pysam

from woldrnaseq.common import (
    add_metadata_arguments,
    add_debug_arguments,
    configure_logging
)
from woldrnaseq.models import (
    load_library_tables,
    load_experiments,
    find_library_bam_file
)

LOGGER = logging.getLogger('make_aggregate_bam')


def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('-n', '--experiment-name', required=True,
                        help='Experiment name to select')
    add_metadata_arguments(parser)
    add_debug_arguments(parser)
    args = parser.parse_args(cmdline)

    configure_logging(args)

    header_printed = False
    libraries = load_library_tables(args.libraries)
    experiments = load_experiments(args.experiments)

    replicates = experiments.loc[args.experiment_name, 'replicates']

    for i, (library_id, library) in enumerate(libraries.loc[replicates].iterrows()):
        filename = find_library_bam_file(library)
        LOGGER.info('  Reading %s %d/%d', filename, i+1, len(replicates))

        mode = get_mode(filename, 'r')
        with pysam.AlignmentFile(filename, mode) as alignment:
            if not header_printed:
                print(str(alignment.header))
                header_printed = True

            for read in alignment:
                print(read.to_string())


def get_mode(alignment_file, mode):
    mode = list(mode)
    if alignment_file.endswith('.sam'):
        mode = mode[:1]

    elif alignment_file.endswith('.bam'):
        if len(mode) == 1:
            mode.append('b')
        else:
            mode[1] = 'b'
    else:
        _, extension = os.path.splitext(alignment_file)
        raise ValueError('Unrecognized filename extension: {}'.format(extension))

    return ''.join(mode)


if __name__ == '__main__':
    main()
