#!/usr/bin/python3
from argparse import ArgumentParser
import os
import time
import re
from glob import glob
import pysam

tenx_root = '/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/'

tenx_bam = {
    1: os.path.join(tenx_root, 'Wold10x-1-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    3: os.path.join(tenx_root, 'Wold10x-3-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    4: os.path.join(tenx_root, 'Wold10x-4-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    5: os.path.join(tenx_root, 'Wold10x-5-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    6: os.path.join(tenx_root, 'Wold10x-6-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    7: os.path.join(tenx_root, 'Wold10x-7-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    8: os.path.join(tenx_root, 'Wold10x-8-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    9: os.path.join(tenx_root, 'Wold10x-9-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    10: os.path.join(tenx_root, 'Wold10x-10-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    11: os.path.join(tenx_root, 'Wold10x-11-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    12: os.path.join(tenx_root, 'Wold10x-12-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    13: os.path.join(tenx_root, 'Wold10x-13-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
}

for cluster in tenx_bam:
    if not os.path.exists(tenx_bam[cluster]):
        print('missing {}'.format(tenx_bam[cluster]))


def get_cluster_barcodes(cluster):
    results = []
    run_re = re.compile(r'barcodes-10x-cluster-{cluster}-run-(?P<run>[\d]+)\.txt'.format(cluster=cluster))
    pattern = 'barcodes-10x-cluster-{cluster}-run-*.txt'.format(cluster=cluster)
    for filename in glob(pattern):
        match = run_re.match(filename)
        run = int(match.group('run'))
        results.append((run, filename))
    return sorted(results)


def get_total_reads(possorted_bam):
    total = 0
    for stat in possorted_bam.get_index_statistics():
        total += stat.total
    return total


def build_barcode_re(barcode_filename):
    barcodes = []
    with open(barcode_filename, 'rt') as instream:
        for line in instream:
            barcode = line.rstrip().split(':')[2]
            barcodes.append('(' + barcode + ')')
    return re.compile('|'.join(barcodes))


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    output_bam = None
    for run, barcode_file in get_cluster_barcodes(args.cluster):
        print('Processing run {} for cluster {}'.format(run, args.cluster))
        barcode_re = build_barcode_re(barcode_file)

        possorted_bam = pysam.AlignmentFile(tenx_bam[run], 'rb')
        if output_bam is None:
            output_bam = pysam.AlignmentFile(args.output, "wb", template=possorted_bam)

        copy_filtered_reads(possorted_bam, output_bam, barcode_re)
        possorted_bam.close()

    if output_bam is not None:
        output_bam.close()


def copy_filtered_reads(possorted_bam, output_bam, barcode_re):
    t0 = time.monotonic()
    tn = t0
    last_percent = 0
    total_reads = get_total_reads(possorted_bam)
    seen = set()
    current_reference = None

    for i, read in enumerate(possorted_bam):
        if read.has_tag('CB'):
            cb = read.get_tag('CB')
            if barcode_re.match(cb) and read.has_tag('UB'):
                ub = read.get_tag('UB')
                unique_key = (cb, ub, read.reference_name, read.reference_start)
                # reset seen dictionary for each new chromosome/reference sequence
                if read.reference_name != current_reference:
                    seen = set()
                    current_reference = read.reference_name
                if unique_key not in seen:
                    output_bam.write(read)
                    seen.add(unique_key)

        next_percent = int(i / total_reads * 100)
        if next_percent > last_percent:
            tn1 = time.monotonic()
            print(next_percent, tn1-tn)
            tn = tn1
            last_percent = next_percent


def make_parser():
    parser = ArgumentParser()
    parser.add_argument('-c', '--cluster', help='cluster to process')
    parser.add_argument('-o', '--output', help='output bam file', required=True)
    return parser


if __name__ == '__main__':
    main()
