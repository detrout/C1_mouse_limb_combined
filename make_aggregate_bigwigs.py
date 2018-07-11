#!/usr/bin/python3
from argparse import ArgumentParser
import os
from pprint import pprint
from collections import Counter, OrderedDict
import logging
import pysam
import pyBigWig

LOGGER = logging.getLogger('make_aggregate_bigwigs')

def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('alignment_files', nargs='+')
    parser.add_argument('--reference-prefix', default='chr')
    parser.add_argument('-o', '--output', required=True,
                        help='Output bigwig filename')
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('-vv', '--debug', action='store_true', default=False)
    
    args = parser.parse_args(cmdline)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARN)
    

    bigwig = pyBigWig.open(args.output, 'w')
    chrom_info = get_chrom_info(args.alignment_files[0], args.reference_prefix)
    bigwig.addHeader(list(chrom_info.items()), maxZooms=10)

    for contig in chrom_info:
        LOGGER.info('Processing: %s', contig)
        current = sum_alignments(args.alignment_files, chrom_info, contig)
        write_contig_to_bigwig(bigwig, contig, current)

    bigwig.close()

def sum_alignments(alignment_files, chrom_info, contig):
    length = chrom_info[contig]
    current = Counter()

    for filename in alignment_files:
        LOGGER.info('  Reading %s', filename)
        current_chrom_info = get_chrom_info(filename)
        assert current_chrom_info[contig] == chrom_info[contig]

        mode = get_mode(filename, 'r')
        with pysam.AlignmentFile(filename, mode) as alignment:
            for reads_at in alignment.pileup(contig, 0, length):
                current[reads_at.pos] += reads_at.n

    return current

def write_contig_to_bigwig(bigwig, contig, counter):
    LOGGER.info('  building bedgraph for %s', contig)
    previous_i = None
    current_value = None
    starts = []
    ends = []
    values = []
    for i in counter:
        if len(starts) == 0:
            current_value = counter[i]
            starts.append(i)
            values.append(float(current_value))
        elif current_value != counter[i]:
            ends.append(previous_i + 1)
            starts.append(i)
            current_value = counter[i]
            values.append(float(current_value))
            
        previous_i = i
    ends.append(previous_i + 1)
    chromosomes = [contig] * len(starts)
    assert len(chromosomes) == len(starts)
    assert len(starts) == len(ends)
    assert len(ends) == len(values)
    bigwig.addEntries(chromosomes, starts, ends=ends, values=values, validate=False)
    LOGGER.info('  wrote %s entries', len(starts))
        
def get_chrom_info(filename, prefix=None):
    chrom_info = OrderedDict()
    mode = get_mode(filename, 'r')
    with pysam.AlignmentFile(filename, mode) as bam:
        for contig, length in zip(bam.references, bam.lengths):
            if prefix is None or contig.startswith(prefix):
                chrom_info[contig] = length

    return chrom_info


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
        _, extension = os.path.splitext(alignemnt_file)
        raise ValueError('Unrecognized filename extension: {}'.format(extension))

    return ''.join(mode)


if __name__ == '__main__':
    main()
