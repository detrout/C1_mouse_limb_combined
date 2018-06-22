#!/usr/bin/python3
"""Figure out what's actually been submitted
"""
from argparse import ArgumentParser
import itertools
import os
import pandas
import requests

from htsworkflow.submission.encoded import ENCODED
from pandasodf import ODFReader

from woldrnaseq import models

from generate_combined_transcript_C1 import paper_433_experiment_files, ASOF_RUN17_experiment_files

def to_files(file_list):
    return [ os.path.expanduser(x.strip()) for x in file_list.split()]

def parse_replicates(replicates):
    for experiment in replicates:
        for library in experiment:
            yield library.replace('_mm10', '').replace('_clean', '')

def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('-s', '--sheet', default=0, help='Sheet to use')
    parser.add_argument('--header', default=None, help="header row")
    parser.add_argument('filename', nargs=1, help='spreadsheet to look at')
    args = parser.parse_args(cmdline)

    header = int(args.header) if args.header is not None else None
    book = ODFReader(args.filename[0])
    data = book.parse(args.sheet, header=header)

    server = ENCODED('www.encodeproject.org')
    server.load_netrc()

    first_experiments = models.load_experiments(to_files(paper_433_experiment_files))
    all_experiments = models.load_experiments(to_files(ASOF_RUN17_experiment_files))

    first_libraries = set(parse_replicates(first_experiments['replicates']))
    all_libraries = set(parse_replicates(all_experiments['replicates']))

    #print(first_libraries)
    #print(all_libraries)
    results = []
    for i, library_id in enumerate(data[data.columns[0]]):
        if library_id in first_libraries:
            tranche = 1
        elif library_id in all_libraries:
            tranche = 2
        else:
            tranche = 'C'

        row = find_library_info(server, library_id)
        row['tranche'] = tranche
        results.append(row)

        if (i+1) % 10:
            print('.', end='', flush=True)

    df = pandas.DataFrame(results)
    df.to_csv('tranche.csv', index=False)


def find_library_info(server, library_id):
    seed = 'barbara-wold:{}'.format(library_id)

    row = {
        'library_id': library_id,
    }

    try:
        data = server.get_json(seed)
        library_accession = data['accession']

        row['library_aliases'] = data['aliases']

        data = server.search_jsonld(searchTerm=library_accession)
        graph = data['@graph']
        if len(graph) > 1:
            raise RuntimeError(
                'Confused: multiple hits for {} {}'.format(library_id, library_accession))
        elif len(graph) == 1:
            row.update({
                'experiment_id': graph[0]['@id'],
                'experiment_alias': graph[0]['aliases']
            })
    except requests.HTTPError:
        pass

    return row
    
if __name__ == '__main__':
    main()
