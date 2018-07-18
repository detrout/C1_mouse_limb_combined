#!/usr/bin/python3
from argparse import ArgumentParser
import os
import pandas
from generate_combined_transcript_C1 import ASOF_RUN17_experiment_files, ASOF_RUN17_library_files
from woldrnaseq import models

def main(cmdline=None):
    

    experiment_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_experiment_files.split() ]
    library_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_library_files.split() ]

    experiments = models.load_experiments(experiment_files)
    libraries = models.load_library_tables(library_files)

    to_include = read_peng_20180710_cluster_memberships()

    found = 0
    missing = 0
    clustered_wigs = {}
    for suffix in ["-mm10-M4-male_all.bw", "-mm10-M4-male_uniq.bw"]:
        for name, cluster in to_include:
            analysis_dir = libraries.loc[name]['analysis_dir']

            wigdir = get_wigdir(analysis_dir)
            wigname = os.path.join(wigdir, name + suffix)
            wigfound = os.path.exists(wigname)
            if wigfound:
                found += 1
                clustered_wigs.setdefault(cluster, []).append(wigname)
            else:
                missing += 1
                print(name, cluster, wigname, wigfound)

        for cluster in clustered_wigs:
            bigwig_name = 'C1_peng_20180710_cluster_bigwigs/' + cluster.replace(' ', '_') + suffix
            args = ['python3', 'merge_bw.py', '-o', os.path.abspath(bigwig_name)]
            args.extend(clustered_wigs[cluster])
            print(' '.join(args))


def get_wigdir(analysis_dir):
    sau_root = '/woldlab/castor/home/sau/flowcells/'
    diane_root = '/woldlab/loxcyc/home/diane/proj/'
    if analysis_dir.startswith(sau_root):
        base, subdir = os.path.split(analysis_dir)
        base = base.replace('C1_e10.5_mouse_limb_run2_June20_2016', 'C1-single-072616')
        return base.replace(sau_root, '/woldlab/castor/home/sau/public_html/')
    elif analysis_dir.startswith(diane_root):
        return analysis_dir.replace(diane_root, '/woldlab/loxcyc/home/diane/public_html/')
    return analysis_dir


def read_peng_20180710_cluster_memberships():
    cluster_membership = []
    with open('C1_peng_20180710_cluster_memberships.txt', 'rt') as instream:
        for line in instream:
            quotes = []
            start = 0
            while True:
                q = line.find("'", start)
                if q == -1:
                    break
                else:
                    start = q + 1
                    quotes.append(q)
            assert len(quotes) == 4
            cell_id = line[quotes[0]+1:quotes[1]]
            if cell_id[0] == 'x':
                cell_id = cell_id[1:]
            value = line[quotes[2]+1:quotes[3]]
            close_brace = line.find("}", quotes[3])
            rest = line[close_brace+1:].strip()
            rest_values = rest.split()
            color = tuple([float(x) for x in rest_values[:3]])
            cluster_name = ' '.join(rest_values[3:])
            cluster_membership.append({
                'cell_id': cell_id,
                'value': value,
                'color': color,
                'cluster_name': cluster_name,
            })
        return cluster_membership

if __name__ == '__main__':
    main()
