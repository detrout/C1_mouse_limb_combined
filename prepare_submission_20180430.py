#!/usr/bin/python3
from argparse import ArgumentParser

import os
import collections
from lxml.html import fromstring
import json
import re
import requests
import glob
import pandas
from urllib.parse import urljoin

from rdflib import Graph, Literal, URIRef

from generate_combined_transcript_C1 import (
    paper_433_experiment_files,
    ASOF_RUN17_experiment_files
)

from woldrnaseq.models import load_experiments

from htsworkflow.util.opener import autoopen
from htsworkflow.util.rdfns import (
    libraryOntology,
    RDF,
    RDFS,
)
from htsworkflow.util.rdfhelp import (
     dump_model,
)

# 20031-20038 are good on flowcell HF7NTBCX2
# 20026-20030 are mixed on flowcell HF7NTBCX2


def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('--first-tranche', default=False, action='store_true',
                        help='Use just the first tranche as experiment list')
    parser.add_argument('--name', required=True, help='submission name')
    parser.add_argument('-s', '--sheet', default=0, help='Sheet to use')
    parser.add_argument('--header', default=None, help="header row")
    parser.add_argument('filename', nargs=1, help='driver spreadsheet')
    args = parser.parse_args(cmdline)
    root_fastq_url = 'http://jumpgate.caltech.edu/runfolders/volvox02/'
    desplit = os.path.expanduser('~/proj/htsworkflow/htsworkflow/pipelines/desplit_fastq.py')

    header = int(args.header) if args.header is not None else None
    data = read_spreadsheet(args.filename[0], args.sheet, header)
    print(data.shape)

    if args.first_tranche:
        experiment_file_list = paper_433_experiment_files.split('\n')
    else:
        experiment_file_list = ASOF_RUN17_experiment_files.split('\n')
    experiment_files = [ os.path.expanduser(x.strip()) for x in  experiment_file_list]
    experiments = load_experiments(experiment_files)
    experiments['replicates'] = experiments['replicates'].apply(lambda l: [x.replace('_mm10', '').replace('_clean', '') for x in l])

    current_experiments = find_experiments_to_submit(experiments, data)

    aliases_tsv = '{}-aliases.tsv'.format(args.name)
    make_library_aliases(current_experiments, aliases_tsv)

    submission_fastqs_tsv = '{}-fastqs.tsv'.format(args.name)
    if not os.path.exists(submission_fastqs_tsv):
        fastq_urls = find_all_fastqs(root_fastq_url, current_experiments, submission_fastqs_tsv)

    fastq_urls = pandas.read_csv(submission_fastqs_tsv, sep='\t')

    barcodes_tsv = '{}-barcodes.tsv'.format(args.name)
    make_library_barcodes(fastq_urls, barcodes_tsv)

    metadata_tsv = '{}-flowcell-details.tsv'.format(args.name)
    metadata = make_metadata(fastq_urls, root_fastq_url, metadata_tsv)

    merge_file = '{}-merge-fastqs.condor'.format(args.name)
    make_desplit_condor(fastq_urls, metadata, desplit, root_fastq_url, merge_file)


def read_spreadsheet(filename, sheet, header=None):
    if filename.endswith('xlsx'):
        data = pandas.read_excel(
            'Second_set_of_limb_single_cell_data_for_Diane_almost_complete_April13_2018.xlsx',
            sheet=sheet,
            header=header
        )
    elif filename.endswith('ods'):
        from pandasodf import ODFReader

        book = ODFReader(filename)
        data = book.parse(sheet, header=header)

    return data


def find_all_fastqs(root_fastq_url, experiments, output_file):
    """Get urls to the raw fastq files for all our replicates
    """
    runfolder = Runfolder(root_fastq_url)

    records = []
    multi = []
    for record in find_replicate_flowcells(experiments):
        fastqs = []
        for flowcell in record['flowcells']:
            fastqs.extend(list(runfolder.find_fastqs(flowcell, record['library_id'])))
        record['fastq_urls'] = fastqs
        fluidigm_fields = parse_fluidigm(urljoin(root_fastq_url, fastqs[0]))
        record['barcode'] = fluidigm_fields['barcode']
        record['location'] = fluidigm_fields['location']
        records.append(record)
        if len(record['flowcells']) > 1:
            multi.append(record)

    df = pandas.DataFrame(records)
    df.to_csv(output_file, sep='\t', index=False)
    if len(multi) > 0:
        print('Warning, runs on multiple flowcells check multiple_flowcells.tsv')
        pandas.DataFrame(multi).to_csv('multiple_flowcells.tsv', sep='\t')

    return df

def find_replicate_flowcells(experiments):
    model = Graph()

    for i, row in experiments.iterrows():
        for extended_id in row.replicates:
            library_id, location, *_ = extended_id.split('_')
            extended_id = library_id + '_' + location
            uri = URIRef('https://felcat.caltech.edu/library/{}/'.format(library_id))
            s = (uri, RDF['type'], libraryOntology['Library'])
            if s not in model:
                model.parse(source=uri, format='rdfa')

            flowcells = model.query("""PREFIX libns: <http://jumpgate.caltech.edu/wiki/LibraryOntology#>

select distinct ?flowcell_id
where {
  ?library a libns:Library ;
           libns:has_lane ?lane .
  ?lane libns:flowcell ?flowcell .
  ?flowcell libns:flowcell_id ?flowcell_id .

}
            """, initBindings={'library': uri})

            yield {'experiment': row.name,
                   'library_id': extended_id,
                   'flowcells': sorted([x[0].value for x in flowcells])
                   }

def find_experiments_to_submit(experiments, submission_table):
    to_upload = set(submission_table[submission_table.columns[0]])
    missing = set(submission_table[submission_table.columns[0]])

    tosubmit = []
    for i, row in experiments.iterrows():
        current = to_upload.intersection(set(row.replicates))
        missing = missing.difference(set(row.replicates))
        if len(current) > 0:
            tosubmit.append({
                'name': row.name,
                'analysis_dir': row.analysis_dir,
                'replicates': list(current)
            })

    print('Not found:', len(missing), sorted(missing))
    df = pandas.DataFrame(tosubmit)
    df.set_index('name', inplace=True)
    return df


def find_seans_fastqs(experiments):
    for i, row in experiments.iterrows():
        for library_id in row.replicates:
            pattern = os.path.join(row.analysis_dir, library_id + '*.fastq.gz')
            files = glob.glob(pattern)
            assert len(files) > 0
            filesets.setdefault(i, []).extend(files)

    #make_desplit_condor(filesets)


def make_library_aliases(experiments, aliases_tsv):
    aliases = {}
    for i, row in experiments.iterrows():
        for library_id in row.replicates:
            aliases.setdefault(row.name, []).append('barbara-wold:{}'.format(library_id))

    with open(aliases_tsv, 'wt') as outstream:
        for key in sorted(aliases):
            outstream.write(key)
            outstream.write('\t')
            outstream.write(','.join(sorted(aliases[key])))
            outstream.write(os.linesep)


def make_library_barcodes(experiments, barcode_tsv):
    def sorted_plate_key(row):
        return row['plate_id'] + '_' + row['plate_location']

    barcodes = {}
    for i, row in experiments.iterrows():
        plate_id, location, *_ = row.library_id.split('_')
        record = {'barcode': row.barcode, 'plate_id': plate_id, 'plate_location': location}
        barcodes.setdefault(row.experiment, []).append(record)

    with open(barcode_tsv, 'wt') as outstream:
        for key in sorted(barcodes):
            outstream.write(key)
            outstream.write('\t')
            outstream.write(json.dumps(sorted(barcodes[key], key=sorted_plate_key)))
            outstream.write(os.linesep)


def make_desplit_condor(experiments, metadata, desplit_cmd, root_url, condor_file):
    """Make condor file to build merged fastqs

    :Parameters:
      - experiments: (pandas.DataFrame) Experiments and their fastq urls
        from find_all_fastqs()
      - metadata: (pandas.DataFrame) metadata details about each fastq
      - desplit_cmd: (filename) Path to the desplit_fastq.py file from htsworkflow
      - condor_file: (filename) target to write condor file
    :Returns:
      True if all the merged fastqs exists, otherwise False
    """
    header = """universe=vanilla
executable=/usr/bin/python3
error=log/desplit_fastq.$(process).out
output=log/desplit_fastq.$(process).out
log=log/desplit_fastq.log
environment="PYTHONPATH=/woldlab/loxcyc/home/diane/proj/htsworkflow"

"""

    experiment_fastqs = {}
    for i, row in metadata.iterrows():
        rl_name = row.experiment + '_' + str(row.read_length)
        output_name = rl_name + '.fastq.gz'
        experiment_fastqs.setdefault(output_name, []).append(row.fastq_url)

    # chunk all fastqs by experiment
    body = []
    for output_name in experiment_fastqs:
        print(output_name)
        fastq_urls = experiment_fastqs[output_name]
        body.extend(['arguments="{} --gzip -o {} {}"'.format(desplit_cmd, output_name, ' '.join(sorted(fastq_urls))),
                     'queue',
                     ''])

    if len(body) > 0:
        with open(condor_file, 'wt') as outstream:
            outstream.write(header)
            outstream.write(os.linesep.join(body))
        return False
    else:
        return True


def make_metadata(experiments, root_fastq_url, filename):
    model = Graph()

    metadata = []
    for i, row in experiments.iterrows():
        fastq_urls = [ urljoin(root_fastq_url, x[1:-1]) for x in row.fastq_urls[1:-1].split(', ')]
        for fastq_url in fastq_urls:
            fastq_data = parse_fluidigm(fastq_url)
            metadata.append({
                'experiment': row.experiment,
                'fastq_url': fastq_url,
                'machine': 'http://jumpgate.caltech.edu/sequencer/8',
                'flowcell': fastq_data['flowcell_id'],
                'lane': fastq_data['lane_number'],
                'barcode': fastq_data['barcode'],
                'read_length': fastq_data['read_length']
            })

    metadata = sorted(metadata, key=lambda row: (row['experiment'], row['flowcell'], row['barcode']))
    df = pandas.DataFrame(metadata, columns=['experiment', 'fastq_url', 'machine', 'flowcell', 'lane', 'barcode', 'read_length'])
    print(df.head())
    df.to_csv(filename, sep='\t', index=False)
    return df


fluidigm_fields = ['library_id', 'location', 'barcode', 'lane_number', 'read']

def parse_fluidigm(pathname):
    path, name = os.path.split(pathname)
    p = r'(?P<library_id>[0-9]{5})_'\
        '(?P<location>[A-H][0-9]{1,2})_'\
        '(?P<barcode>[AGCT-]+)_'\
        'L00(?P<lane_number>[1-8])_'\
        'R(?P<read>[1-3])'

    match = re.match(p, name)
    if match is not None:
        fields = { k: match.group(k) for k in fluidigm_fields }
        with autoopen(pathname, 'rt') as stream:
            fields.update(parse_fastq_header(stream.readline()))
            seq = stream.readline()
            fields['read_length'] = len(seq)
        return fields

def parse_fastq_header(header):
    header = header.strip()
    read_id, extra = header.split(' ')
    fields = read_id.split(':')
    extra_fields = extra.split(':')
    return {
        'flowcell_id': fields[2],
        #'lane_number': fields[3],
        #'read': fields[4],
        'barcode': extra_fields[3],
    }

class Runfolder:
    def __init__(self, root_url):
        self.root_url = root_url
        self.pages = {}

    def load_index(self, url=''):
        absolute_url = urljoin(self.root_url, url)

        response = requests.get(absolute_url)
        if response.status_code != 200:
            raise RuntimeError('Unable to access {}. Status {}'.format(absolute_url, response.status_code))
        tree = fromstring(response.content)
        rows = tree.xpath('*/table/tr/td/a')
        if len(rows) == 0:
            raise RuntimeError('{} is not a directory'.format(absolute_url))
        if rows[0].text == 'Parent Directory':
            rows.pop(0)
        self.pages[url] = [ x.text for x in rows ]

    def find_flowcell(self, flowcell):
        root = ''
        if root not in self.pages:
            self.load_index(root)

        for name in self.pages[root]:
            if flowcell in name:
                return name

    def _find_unaligned(self, url):
        if url not in self.pages:
            self.load_index(url)

        for name in self.pages[url]:
            for unaligned in ['Unaligned.dualIndex/', 'Unaligned/']:
                if unaligned == name:
                    return url + name

        raise RuntimeError('Unable to find index in {}'.format(url))

    def _find_extended_id(self, url, extended_id):
        if url not in self.pages:
            self.load_index(url)

        for name in self.pages[url]:
            if extended_id in name:
                return url + name

    def find_fastqs(self, flowcell, extended_id):
        runfolder = self.find_flowcell(flowcell)
        assert runfolder is not None
        unaligned = self._find_unaligned(runfolder)
        assert unaligned is not None
        project = self._find_extended_id(unaligned, extended_id)
        sample = self._find_extended_id(project, extended_id)

        if sample not in self.pages:
            self.load_index(sample)

        for name in self.pages[sample]:
            if 'fastq.gz' in name:
                yield sample + name

if __name__ == '__main__':
    main()
