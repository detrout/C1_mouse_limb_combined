#!/usr/bin/python3
"""Create publication set for the C1 mouse paper
"""
import logging
import pandas
from pprint import pprint
from htsworkflow.submission.encoded import ENCODED, DCCValidator

def main():
    server = ENCODED('www.encodeproject.org')
    server.load_netrc()
    validator = DCCValidator(server)

    table = pandas.read_csv('publication_files.tsv', sep='\t', skipfooter=1)

    award = '/awards/UM1HG009443/'
    lab = '/labs/barbara-wold/'

    payload = {
        'award': award,
        'lab': lab,
        'related_files': list(table['file']),
        'references': ['/publications/e0d01543-9965-4edb-933c-778a40575cd9/'],
    }

    #validator.validate(payload, 'publication_data')
    print(server.post_json('/publication-data/', payload))
    #pprint(payload)

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main()
