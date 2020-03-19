#!/usr/bin/python3
from argparse import ArgumentParser
import csv
from collections import namedtuple
import os
import pysam

tenx_root = '/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/'

tenx_bam = {
    '1': os.path.join(tenx_root, 'Wold10x-1-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '3': os.path.join(tenx_root, 'Wold10x-3-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '4': os.path.join(tenx_root, 'Wold10x-4-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '5': os.path.join(tenx_root, 'Wold10x-5-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '6': os.path.join(tenx_root, 'Wold10x-6-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '7': os.path.join(tenx_root, 'Wold10x-7-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '8': os.path.join(tenx_root, 'Wold10x-8-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '9': os.path.join(tenx_root, 'Wold10x-9-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '10': os.path.join(tenx_root, 'Wold10x-10-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '11': os.path.join(tenx_root, 'Wold10x-11-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '12': os.path.join(tenx_root, 'Wold10x-12-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
    '13': os.path.join(tenx_root, 'Wold10x-13-encode-count-cells10000', 'outs', 'possorted_genome_bam.bam'),
}

for cluster in tenx_bam:
    if not os.path.exists(tenx_bam[cluster]):
        print('missing {}'.format(tenx_bam[cluster]))


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    clusters = build_cluster_recognizer('../monocle/mouse/barcodes-to-cluster.csv')
    possorted_bam = pysam.AlignmentFile(tenx_bam[args.run], 'rb')
    if args.header:
        print(str(possorted_bam.header))

    filter_reads(possorted_bam, clusters[args.cluster][args.run])
    possorted_bam.close()


def make_parser():
    parser = ArgumentParser()
    parser.add_argument('-c', '--cluster', help='cluster to process')
    parser.add_argument('-r', '--run', help='run number to process')
    parser.add_argument('--header', help='print bam file header', default=False, action='store_true')
    return parser


def filter_reads(possorted_bam, barcodes):
    seen = set()
    current_reference = None

    for i, read in enumerate(possorted_bam):
        if read.has_tag('CB'):
            cb = read.get_tag('CB')
            if is_recognized_barcode(barcodes, cb) and read.has_tag('UB'):
                ub = read.get_tag('UB')
                unique_key = (cb, ub, read.reference_name, read.reference_start)
                # reset seen dictionary for each new chromosome/reference sequence
                if read.reference_name != current_reference:
                    seen = set()
                    current_reference = read.reference_name
                if unique_key not in seen:
                    print(read.to_string())
                    seen.add(unique_key)


CellBarcode = namedtuple('CellBarcode', ['run', 'stage', 'barcode', 'cluster'])


def parse_barcode_csv(filename, limit=None):
    with open(filename) as barcodes:
        reader = csv.reader(barcodes)
        next(reader)
        for i, row in enumerate(reader):
            run = []
            stage = []
            field0 = row[0]
            run_ = field0.find('_') + 1
            for c in field0[4:run_-1]:
                run.append(c)
            # +2 to be one after the trailing 0 or 5
            run_2 = field0.find('_', run_) + 2
            for c in field0[run_:run_2]:
                if c == '_':
                    c = '.'
                stage.append(c)
            run = ''.join(run)
            stage = ''.join(stage)

            barcode = field0[run_2:-2] + '-' + field0[-1:]
            cluster = row[1]
            yield CellBarcode(run, stage, barcode, cluster)
            if limit is not None and i > limit:
                break


def build_cluster_recognizer(filename):
    clusters = {}

    for row in parse_barcode_csv(filename):
        run = clusters.setdefault(row.cluster, {})
        barcode_tree = run.setdefault(row.run, {})
        barcode_node = barcode_tree
        for c in row.barcode:
            barcode_node = barcode_node.setdefault(c, {})

    return clusters


def is_recognized_barcode(barcode_tree, cell_barcode):
    barcode_node = barcode_tree
    for c in cell_barcode:
        try:
            barcode_node = barcode_node[c]
        except KeyError:
            return False
    return True


if __name__ == '__main__':
    main()
