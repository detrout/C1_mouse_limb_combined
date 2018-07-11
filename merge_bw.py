#!/usr/bin/python3
from argparse import ArgumentParser
import os
from subprocess import check_call, check_output, PIPE
from hashlib import md5

def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('alignment_files', nargs='+')
    parser.add_argument('--reference-prefix', default='chr')
    parser.add_argument('-o', '--output', required=True,
                        help='Output bigwig filename')
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('-vv', '--debug', action='store_true', default=False)
    args = parser.parse_args(cmdline)

    ucsc_tools = '/woldlab/castor/proj/programs/x86_64-366/'
    bigWigMerge = os.path.join(ucsc_tools, 'bigWigMerge')
    bedSort = os.path.join(ucsc_tools, 'bedSort')
    bedGraphToBigWig = os.path.join(ucsc_tools, 'bedGraphToBigWig')
    chrom_sizes = os.path.expanduser('~/proj/C1_mouse_limb_combined/mm10.chrom.sizes')

    uniq_name = md5(b' '.join([ x.encode('ascii') for x in args.alignment_files])).hexdigest()
    merged_bg = uniq_name + '-merged.bg'
    filtered_bg = uniq_name + '-filtered.bg'
    sorted_bg = uniq_name + '-sorted.bg'
    
    check_call([bigWigMerge] + args.alignment_files + [merged_bg])
    check_output('grep chr ' + merged_bg + ' > ' + filtered_bg, shell=True)
    os.unlink(merged_bg)
    check_call([bedSort, filtered_bg, sorted_bg])
    os.unlink(filtered_bg)
    check_call([bedGraphToBigWig, sorted_bg, chrom_sizes, args.output])
    os.unlink(sorted_bg)

if __name__ == '__main__':
    main()
