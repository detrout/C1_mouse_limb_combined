#!/usr/bin/python3
from argparse import ArgumentParser
import os
import pandas
import subprocess

from woldrnaseq.models import load_library_tables, genome_name_from_library


def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('libraries', nargs='+')
    args = parser.parse_args(cmdline)

    results = []
    for library in args.libraries:
        results.append(build_hash_tree(library))

    results = pandas.DataFrame(results)
    results.to_csv('compare_bams_genome.csv')


def build_hash_tree(library_filename):
    table = load_library_tables([library_filename])

    hashes = {}
    for library_id, row in table.iterrows():
        analysis_dir = row.analysis_dir
        name = row.analysis_name + '-' + genome_name_from_library(row) + '_genome.bam'
        alignment = os.path.join(analysis_dir, name)
        hashes[library_id] = hash_alignments(alignment)

    return hashes

def hash_alignments(filename):
    samtools = subprocess.Popen(['samtools', 'view', filename], stdout=subprocess.PIPE)
    sha = subprocess.Popen(['sha256sum'], stdin=samtools.stdout, stdout=subprocess.PIPE)
    samtools.stdout.close()
    output = sha.communicate()[0]
    return output.decode('ascii')[:-4]

if __name__ == '__main__':
    main()
