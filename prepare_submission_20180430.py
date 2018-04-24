import os
import glob
import pandas

from generate_combined_transcript_C1 import ASOF_RUN17_experiment_files

from woldrnaseq.models import load_experiments

def main():

    desplit = os.path.expanduser('~/proj/htsworkflow/htsworkflow/pipelines/desplit_fastq.py')
    data = pandas.read_excel(
        'Second_set_of_limb_single_cell_data_for_Diane_almost_complete_April13_2018.xlsx',
        sheet=0,
        header=None
    )
    experiment_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_experiment_files.split('\n') ]
    experiments = load_experiments(experiment_files)

    to_upload = set(data[0])
    missing = set(data[0])
    filesets = {}
    for i, row in experiments.iterrows():
        current = to_upload.intersection(set(row.replicates))
        missing = missing.difference(set(row.replicates))
        for library_id in current:
            pattern = os.path.join(row.analysis_dir, library_id + '*.fastq.gz')
            files = glob.glob(pattern)
            assert len(files) > 0
            filesets.setdefault(i, []).extend(files)

    #make_desplit_condor(filesets)

def make_desplit_condor(filesets):
    print("""universe=vanilla
executable=/usr/bin/python3
error=log/desplit_fastq.$(process).out
output=log/desplit_fastq.$(process).out
log=log/desplit_fastq.log
environment="PYTHONPATH=/woldlab/loxcyc/home/diane/proj/htsworkflow"

""")
    for key in filesets:
        output_name = key + '.fastq.gz'
        print('arguments="{} --gzip -o {} {}"'.format(desplit, output_name, ' '.join(sorted(filesets[key]))))
        print('queue')
        print('')

if __name__ == '__main__':
    main()
