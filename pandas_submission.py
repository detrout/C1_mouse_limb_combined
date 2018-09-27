"""Start prototyping scanning using pandas instead of sparql
"""
import argparse
import os
import json
import logging
import pandas
import sys
from htsworkflow.util.hashfile import make_md5sum
from htsworkflow.submission.aws_submission import upload_file
from htsworkflow.submission.fastqname import FastqName
from htsworkflow.submission.encoded import DCCValidator, ENCODED, TYPE_TO_COLLECTION
from htsworkflow.util import opener


PANDAS_ODF = os.path.expanduser('~/src/pandasodf')
if PANDAS_ODF not in sys.path:
    sys.path.append(PANDAS_ODF)
from pandasodf import ODFReader

LOGGER = logging.getLogger('pandas_submission')

def main(cmdline=None):
    logging.basicConfig(level=logging.INFO)

    #server = ENCODED('test.encodedcc.org')
    #spreadsheet_name = 'C1-mouse-forlimb-submission-201804-testserver.ods'

    server = ENCODED('www.encodeproject.org')
    spreadsheet_name = 'C1-mouse-forlimb-submission-201804.ods'
    server.load_netrc()

    dry_run=False
    book = ODFReader(spreadsheet_name)
    process_fastqs(server, book, dry_run)

def process_fastqs(server, book, dry_run):
    award = '/awards/UM1HG009443/'
    lab = '/labs/barbara-wold/'

    files = book.parse('File', header=0)

    LOGGER.info('Attaching additional file metadata')
    files['flowcell_details:json'] = files['submitted_file_name'].apply(add_fastq_metadata)
    files['file_size:integer'] = files['submitted_file_name'].apply(add_file_size)
    files['read_length:integer'] = files['submitted_file_name'].apply(add_read_length)
    files['md5sum'] = files['submitted_file_name'].apply(add_md5s)

    LOGGER.info('Validating metadata')
    validator = DCCValidator(server=server)
    validate(server, validator, book, files)

    LOGGER.info('Uploading files')
    upload(server, validator, files, dry_run=dry_run)

def validate(server, validator, book, files):
    # load background attributes
    for sheet_name in ['Biosample', 'Library', 'Experiment', 'Biosample']:
        collection = TYPE_TO_COLLECTION[sheet_name]
        sheet = book.parse(sheet_name, header=0)
        server.prepare_objects_from_sheet(collection, sheet, validator)

    # validate our annotated files object
    collection = TYPE_TO_COLLECTION['File']
    server.prepare_objects_from_sheet(collection, files, validator)
    print("Validation passed")


def upload(server, validator, files, dry_run=True, retry=False):
    to_create = server.prepare_objects_from_sheet('/files/', files, validator=validator)
    for i, new_object in to_create:
        upload_file(server, validator, new_object, dry_run, retry)


def shorten_filename(submission_pathname):
    _, pathname = os.path.split(submission_pathname)
    return pathname


def make_aliases(short_name):
    return ['barbara-wold:'+short_name]


def add_md5s(submission_pathname):
    LOGGER.debug("Updating file md5sum %s", submission_pathname)
    md5 = make_md5sum(submission_pathname)
    if md5 is None:
        errmsg = "Unable to produce md5sum for {0}"
        LOGGER.warning(errmsg.format(submission_pathname))
    else:
        return md5

def add_file_size(submission_pathname):
    file_size = os.stat(submission_pathname).st_size
    return file_size

def add_read_length(submission_pathname):
    stream = opener.autoopen(submission_pathname, 'rt')
    header = stream.readline().strip()
    sequence = stream.readline().strip()
    stream.close()
    read_length = len(sequence)
    LOGGER.debug("Updating read length: %d", read_length)
    return read_length

def add_fastq_metadata(submission_pathname):
    _, filename = os.path.split(submission_pathname)
    experiment_name = filename.replace('.fastq.gz', '')

    df = pandas.read_csv('submission-201804-flowcell-details.tsv', sep='\t', index_col=False)
    details = df[df['experiment'] == experiment_name]

    cell = []
    for i, row in details.iterrows():
        cell.append({
            'machine': row.machine,
            'flowcell': row.flowcell,
            'lane': str(row.lane),
            'barcode': row.barcode,
        })

    return json.dumps(cell)


if __name__ == '__main__':
    main()
