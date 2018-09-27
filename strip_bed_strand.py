#!/usr/bin/python3
from argparse import ArgumentParser
import os


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    with open(args.infile, 'rt') as instream:
        with open(args.outfile, 'wt') as outstream:
            for line in instream:
                fields = line.rstrip().split('\t')
                fields[5] = '.'
                outstream.write('\t'.join(fields))
                outstream.write(os.linesep)

def make_parser():
    parser = ArgumentParser()
    parser.add_argument('-i', '--infile', help='input bed file')
    parser.add_argument('-o', '--outfile', help='output bed file')
    return parser

if __name__ == '__main__':
    main()
