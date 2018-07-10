#!/usr/bin/python3
from argparse import ArgumentParser
import os
from pprint import pprint
from collections import Counter
import pysam
import pyBigWig

def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('alignment_files', nargs='+')
    parser.add_argument('--reference-prefix', default='chr')
    args = parser.parse_args(cmdline)


    chrom_info = get_chrom_info(args.alignment_files[0], args.reference_prefix)
    for contig in chrom_info:

        current = sum_alignments(args.alignment_files, chrom_info, contig)
        print('current', len(current))
        for i, key  in enumerate(current):
            print(key, current[key])
        break


def sum_alignments(alignment_files, chrom_info, contig):
    length = chrom_info[contig]
    current = Counter()

    for filename in alignment_files:
        current_chrom_info = get_chrom_info(filename)
        assert current_chrom_info[contig] == chrom_info[contig]

        mode = get_mode(filename, 'r')
        with pysam.AlignmentFile(filename, mode) as alignment:
            for reads_at in alignment.pileup(contig, 0, length):
                current[reads_at.pos] += reads_at.n
                #read = next(reads_at)
                #print(read.reference_start)

    return current


def get_chrom_info(filename, prefix=None):
    chrom_info = {}
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
