import logging
import pandas
from pprint import pprint
from htsworkflow.submission.encoded import ENCODED, DCCValidator


def main():
    server = ENCODED('www.encodeproject.org')
    server.load_netrc()

    released_experiments = [
        '/experiments/ENCSR062KGY/',
        '/experiments/ENCSR723FBU/',
        '/experiments/ENCSR559CDN/',
        '/experiments/ENCSR160DGP/',
        '/experiments/ENCSR881ZYX/',
        '/experiments/ENCSR220RKA/',
        '/experiments/ENCSR182LFI/',
        '/experiments/ENCSR716HQB/',
        '/experiments/ENCSR652JLT/',
        '/experiments/ENCSR182WHH/',
        '/experiments/ENCSR839DYB/',
        '/experiments/ENCSR311IKT/',
        '/experiments/ENCSR938RJZ/',
        '/experiments/ENCSR546KIB/',
        '/experiments/ENCSR430OIC/',
    ]


    files = []
    for experiment in released_experiments:
        data = server.get_json(experiment)
        for f in data['files']:
            #print(experiment, f['@id'])
            files.append(f['@id'])

    accession = '/publication-data/ENCSR226XLF/'
    data = server.get_json(accession)
    if len(data['files']) > 50:
        with open('ENCSR226XLF.json', 'wt') as outstream:
            outstream.write('ENCSR226XLF_old.json')

    print(server.patch_json(accession, {'related_files': files}))
    
    data = server.get_json(accession)
    print(data['files'])

if __name__ == '__main__':
    main()
