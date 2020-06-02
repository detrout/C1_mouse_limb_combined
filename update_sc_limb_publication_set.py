#!/usr/bin/python3
"""Update the single cell forelimb publication set

Update the publication set with what we actually submitted.
"""
from htsworkflow.submission.encoded import ENCODED
import pandas


def main(cmdline=None):
    data_id = 'ENCSR226XLF'

    server = ENCODED('www.encodeproject.org')
    server.load_netrc()

    book_name = 'C1-mouse-forelimb-submission-201907-uploaded-production.xlsx'
    book = pandas.ExcelFile(book_name)
    file_sheet = book.parse('File')
    files = ['/files/{}/'.format(x) for x in file_sheet['accession']]
    print('Have {} files to update with'.format(len(files)))

    publication = server.get_json(data_id)
    print(publication['description'])
    print('Publication has {} files'.format(len(publication['files'])))
    p = input("Replace with {} files from {}? ".format(len(files), book_name))

    if p.lower().startswith('y'):
        response = server.patch_json(data_id, {'related_files': files})
        print(response)


if __name__ == "__main__":
    main()

