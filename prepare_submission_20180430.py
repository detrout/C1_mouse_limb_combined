#!/usr/bin/python3
import os
import collections
import re
import glob
import pandas
from urllib.parse import urljoin

from rdflib import Graph, Literal, URIRef

from generate_combined_transcript_C1 import ASOF_RUN17_experiment_files

from woldrnaseq.models import load_experiments

from htsworkflow.util.opener import autoopen
from htsworkflow.util.rdfns import (
    libraryOntology,
)
from htsworkflow.util.rdfhelp import (
     dump_model,
)

# 20031-20038 are good on flowcell HF7NTBCX2
# 20026-20030 are mixed on flowcell HF7NTBCX2

def main():
    root_fastq_url = 'http://jumpgate.caltech.edu/runfolders/volvox02/'
    desplit = os.path.expanduser('~/proj/htsworkflow/htsworkflow/pipelines/desplit_fastq.py')
    data = pandas.read_excel(
        'Second_set_of_limb_single_cell_data_for_Diane_almost_complete_April13_2018.xlsx',
        sheet=0,
        header=None
    )
    experiment_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_experiment_files.split('\n') ]
    experiments = load_experiments(experiment_files)

    current_experiments = find_experiments_to_submit(experiments, data)

    aliases_tsv = 'submission-201804-aliases.tsv'
    make_library_aliases(current_experiments, aliases_tsv)

    #make_desplit_condor(fastq_urls, desplit, root_fastq_url, 'merge_20180430_fastqs.condor')

def find_experiments_to_submit(experiments, submission_table):
    to_upload = set(submission_table[0])
    missing = set(submission_table[0])

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
def make_desplit_condor(experiments, desplit_cmd, root_url, condor_file):
    """Make condor file to build merged fastqs

    :Parameters:
      - experiments: (pandas.DataFrame) Experiments and their fastq urls
        from find_all_fastqs()
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
    for i, row in experiments.iterrows():
        output_name = row.experiment + '.fastq.gz'
        if not os.path.exists(output_name):
            fastq_urls = [ urljoin(root_url, x[1:-1]) for x in row.fastq_urls[1:-1].split(', ')]
            experiment_fastqs.setdefault(output_name, []).extend(fastq_urls)

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


if __name__ == '__main__':
    main()
