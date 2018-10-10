from htsworkflow.submission.encoded import ENCODED
import pandas
from pprint import pprint
from pandasodf import ODFReader

def main():
    data_id = 'ENCSR574CRQ'
    #test_id = 'TSTSR910688'

    server = ENCODED('www.encodeproject.org')
    #server = ENCODED('test.encodedcc.org')
    server.load_netrc()

    experiments = {}
    for e in load_experiments_from_ods('C1-encode3-limb-tranche1-resubmit.ods'):
        experiments['/experiments/{}/'.format(e)] = None
    for e in load_experiments_from_ods('C1-mouse-forlimb-submission-201804.ods'):
        experiments['/experiments/{}/'.format(e)] = None

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

    # load bulk list
    for l in df['Library number']:
        graph = server.search_jsonld(searchTerm='barbara-wold:{}'.format(l))
        for result in graph['@graph']:
            if 'Experiment' in result['@type']:
                experiment = server.get_json(result['@id'])
                if experiment['status'] == 'released':
                    experiments[result['@id']] = None

    pub_id = '/publications/e0d01543-9965-4edb-933c-778a40575cd9/'
    pub = server.get_json(pub_id)
    dataset = pub.get('datasets', [])
    dataset = [ x['@id'] for x in dataset ]
    #for d in dataset:
    #    if d in experiments:
    #        del experiments[d]
    experiments = list(experiments)
    print('posted', len(dataset))
    print('update', len(experiments))
    pprint(list(experiments))
    print(type(experiments))
    print(type(dataset))
    print(set(experiments).difference(set(dataset)))
    print(server.patch_json(pub_id, {'datasets': experiments}))


def load_experiments_from_ods(filename):
    book = ODFReader(filename)
    sheet = book.parse('Experiment')
    return sheet['experiment_accession'].values

if __name__ == '__main__':
    main()
