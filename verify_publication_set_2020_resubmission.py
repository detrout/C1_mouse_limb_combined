#!/usr/bin/python3
"""Verify the publication set after we split everything

I need to update my publication set after I replaced
a bunch of fastqs and experiments when we split every cell
into its own experiment.
"""
from htsworkflow.submission.encoded import ENCODED
import pandas
from pprint import pprint


def main(cmdline=None):
    data_id = 'ENCSR574CRQ'

    server = ENCODED('www.encodeproject.org')
    server.load_netrc()

    book = pandas.ExcelFile('C1-mouse-forelimb-submission-201907-uploaded-production.xlsx')
    files = book.parse('File')

    publication = server.get_json(data_id)
    publication_accessions = [x[len('/files/'):-1] for x in publication['files']]
    print('submitted', len(files['accession']))
    print('posted', len(publication_accessions))
    print('Intersect', len(set(files['accession'].values).intersection(publication_accessions)))


if __name__ == "__main__":
    main()

