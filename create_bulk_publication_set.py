from htsworkflow.submission.encoded import ENCODED
import pandas
from pprint import pprint

def main():
    data_id = 'ENCSR574CRQ'
    #test_id = 'TSTSR910688'

    #libraries = ['17298', '17299', '15019', '15020', '16930',
    #             '16931', '16110', '16111', '15084', '15085',
    #             '16134', '16135',]

    server = ENCODED('www.encodeproject.org')
    #server = ENCODED('test.encodedcc.org')
    server.load_netrc()
    df = pandas.read_excel('Mouse embryo samples list library numbers Diane August 21 2017.xlsx',
                           sheet='Sheet 1',
                           usecols=[0,1,2,3],
                           dtype={
                               2: str,
                               3: str
                           }
    )
    df = df.dropna()
    print(df)

    award = '/awards/UM1HG009443/'
    lab = '/labs/barbara-wold/'

    files = {}
    for l in df['Library number']:
        graph = server.search_jsonld(searchTerm='barbara-wold:{}'.format(l))
        for result in graph['@graph']:
            if 'Experiment' in result['@type']:
                experiment = server.get_json(result['@id'])
                print(result['@id'], l)
                for f in experiment['files']:
                    assembly = f.get('assembly')
                    genome_annotation = f.get('genome_annotation')
                    if assembly is None and genome_annotation is None:
                        files[f['@id']] = f.get('submitted_file_name')
                        print(f['@id'], f.get('submitted_file_name'))
                    elif assembly == 'mm10' and genome_annotation == 'M4':
                        files[f['@id']] = f.get('submitted_file_name')
                        print(f['@id'], f.get('submitted_file_name'))
                    else:
                        print('  ', assembly, genome_annotation)

    payload = {
        'award': award,
        'lab': lab,
        'related_files': list(files.keys()),
        'references': ['/publications/e0d01543-9965-4edb-933c-778a40575cd9/'],
    }
    #pprint(payload)
    #if data_id is None:
    #    result = server.post_json('/publication-data/', payload)
    #    print(result['@id'])
    results = server.patch_json(data_id, {'related_files': list(files)})
    print(results)
    print(results.get('@id'))

if __name__ == '__main__':
    main()
