"""For all libraries, do we have a valiad alias?
"""
import pandas
from requests.exceptions import HTTPError
from htsworkflow.submission.encoded import ENCODED
from to_include import generate_to_include

def main():
    server = ENCODED('www.encodeproject.org')
    server.load_netrc()

    to_include = [x.replace('_mm10', '') for x in generate_to_include()]

    found = []
    for lib in to_include[1:]:
        alias = 'barbara-wold:' + lib
        try:
            obj = server.get_json(alias)
        except HTTPError as e:
            print(lib, 'not found')
        else:
            found.append((lib,obj['@id']))

    df = pandas.DataFrame(found, columns=['library_id', 'id'])
    print(df.head())
    df.to_csv('library_id_to_accession.tsv', sep='\t', index=False)

if __name__ == '__main__':
    main()
