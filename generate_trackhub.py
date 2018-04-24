#!/usr/bin/python3
"""Generate trackhub for specified set...
"""
import argparse
import os
import numpy
import pandas
import trackhub

from woldrnaseq import models
from woldrnaseq.make_tracks import make_bigwig_track_name, make_bam_track_name

from generate_combined_transcript_C1 import ASOF_RUN17_library_files

paper_433_bigwig_root = """~sau/public_html/H5LV3BCXY/
~diane/proj/C1_e10.5_mouse_limb_run3_Dec5_2016_2/
~sau/public_html/C1_mouse_limb_combined_Mar_2017/
~sau/public_html/HVGWNBCXX/
"""

asof_run17_bigwig_paths = paper_433_bigwig_root + """~sau/public_html/HFNLNBCX2/
~sau/public_html/H7CNTBCX2/
~sau/public_html/HFNLTBCX2/
~sau/public_html/HF7NTBCX2/
~sau/public_html/HFNYNBCX2/
"""

cluster_mapping = {
    'deep_red': 'deep_red',
    'red': 'red',
    'black': 'black',
    'green': 'green'
}
# convert human friendly name to machine friendly
human_cluster_mapping = {
    'deep red': 'deep_red',
    'red': 'red',
    'black': 'black',
    'green': 'green'
}
colors = {
    'deep red': '#843c39',
    'red': '#d62728',
    'black': '#000000',
    'green': '#2ca02c',
}
track_dir = 'limb_cells_for_track_hub_splice_isoforms'
public_dir = os.path.expanduser('~/public_html/' + track_dir)

def main(cmdline=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry-run', action='store_true', default=False)

    args = parser.parse_args(cmdline)
    
    combined = 'Limb_cells_for_track_hub_splice_isoforms_red_green_black_only_April17_2018.csv'
    if not os.path.exists(combined):
        roots = [ os.path.expanduser(x.strip()) for x in asof_run17_bigwig_paths.split() ]
        df = load_limb_cells_20180417_set()
        libraries = load_asof_run17_libraries()
        libraries = pandas.merge(df, libraries, how='inner', left_index=True, right_index=True)

        libraries = add_bigwig_paths(libraries, roots)
        libraries = add_bam_paths(libraries)
        libraries.to_csv(combined)
    else:
        libraries = pandas.read_csv(combined)
    
    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name='C1_mouse_limb',
        short_label='C1 mouse limb',
        long_label='C1 mouse limb April 17 2018',
        genome='mm10',
        email='diane@caltech.edu')

    make_bigwig_trackhub(libraries, trackdb)
    make_bam_trackhub(libraries, trackdb)

    if not args.dry_run:
        trackhub.upload.upload_hub(
            hub=hub,
            host='localhost',
            remote_dir=public_dir)
    else:
        print(trackdb)

    #print(hub.hub)
    print('trackhub: ' + 'http://woldlab.caltech.edu/~diane/' + track_dir + '/' + hub.hub + '.hub.txt')
    
def load_limb_cells_20180417_set():
    df = pandas.read_excel(
        'Limb_cells_for_track_hub_splice_isoforms_red_green_black_only_April17_2018.xlsx',
        sheet=0,
    )
    df.columns = ['library_id', 'cluster_assignment']
    df = df.set_index('library_id')
    filtered =  df.dropna(axis=0, how='any')
    print('Filtered to cluster members', filtered.shape)
    return filtered


def load_asof_run17_libraries():
    library_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_library_files.split('\n') ]
    libraries = models.load_library_tables(library_files)
    name = libraries.index.name
    libraries.index = [ x.replace('_mm10', '').replace('_clean', '') for x in libraries.index]
    libraries.index.name = name

    return libraries

def add_bigwig_paths(libraries, roots):
    for track_type in ['uniq', 'all']:
        bigwigs = []
        for library_id, row in libraries.iterrows():
            for root in roots:
                try:
                    track_name = make_bigwig_track_name(row, track_type, root)
                except ValueError:
                    track_name = None
                if track_name is not None:
                    track_name = os.path.join(root, track_name)
                    break
            bigwigs.append(track_name)
        libraries[track_type] = bigwigs
    return libraries

def add_bam_paths(libraries):
    bams = []
    for library_id, row in libraries.iterrows():
        track_name = make_bam_track_name(row, row.analysis_dir)
        track_name = os.path.join(row.analysis_dir, track_name)
        assert os.path.exists(track_name)
            
        bams.append(track_name)
    libraries['bam'] = bams
    return libraries

def make_bigwig_trackhub(libraries, trackdb):

    cluster = trackhub.SubGroupDefinition(
            name='cluster',
            label='cluster',
            mapping=cluster_mapping)

    subgroups = [
        cluster,
        trackhub.SubGroupDefinition(
            name='multi',
            label='multi',
            mapping={
                'uniq': 'uniq',
                'all': 'all',
            }),
    ]

    composite = trackhub.CompositeTrack(
        name='composite',
        short_label='signal',
        dimensions='dimX=cluster dimY=multi',
        tracktype='bigWig',
        visibility='dense',
    )
    composite.add_subgroups(subgroups)
    trackdb.add_tracks(composite)

    signal_view = trackhub.ViewTrack(
        name='signalviewtrack',
        view='signal',
        visibility='dense',
        tracktype='bigWig',
        short_label='Signal')
    composite.add_view(signal_view)

    for i, (library_id, row) in enumerate(libraries.iterrows()):
        for track_index, track_type in enumerate(['uniq', 'all']):
            cluster_assignment = human_cluster_mapping[row.cluster_assignment]
            priority = "{:04d}".format(i*10+track_index)
            track = trackhub.Track(
                url=make_home_url(row[track_type]),
                name=priority + '_' + row.analysis_name + '_' + track_type,
                visibility='dense',
                tracktype='bigWig',
                subgroups={'cluster': cluster_assignment, 'multi': track_type},
                color=hex_to_ucsc_color(colors[row.cluster_assignment]),
                priority=i * 10 + track_index,
            )
            signal_view.add_tracks(track)


def make_bam_trackhub(libraries, trackdb):
    cluster = trackhub.SubGroupDefinition(
            name='cluster',
            label='cluster',
            mapping=cluster_mapping)

    subgroups = [
        cluster,
        trackhub.SubGroupDefinition(
            name='multi',
            label='multi',
            mapping={
                'uniq': 'uniq',
                'all': 'all',
            }),
    ]

    bam_composite = trackhub.CompositeTrack(
        name='reads',
        short_label='reads',
        dimensions='dimX=cluster',
        tracktype='bam',
        visibility='dense',
    )
    bam_composite.add_subgroups([cluster])
    trackdb.add_tracks(bam_composite)
    
    bam_view = trackhub.ViewTrack(
        name='readview',
        view='reads',
        visibility='dense',
        tracktype='bam',
        short_label='Reads')
    bam_composite.add_view(bam_view)

    curdir = os.getcwd()
    os.chdir(public_dir)

    try:
        for i, (library_id, row) in enumerate(libraries.iterrows()):
            cluster_assignment = human_cluster_mapping[row.cluster_assignment]
            priority = "{:03d}".format(i)
            path, name = os.path.split(row.bam)
            bam_index = row.bam + '.bai'
            if not os.path.exists(name):
                os.symlink(row.bam, name)
            if not os.path.exists(name + '.bai'):
                os.symlink(bam_index, name + '.bai')

            url = 'http://woldlab.caltech.edu/~diane/' + track_dir + '/' + name
            track = trackhub.Track(
                url=url,
                name=priority + '_' + row.analysis_name + '_reads',
                visibility='dense',
                tracktype='bam',
                subgroups={'cluster': cluster_assignment},
                color=hex_to_ucsc_color(colors[row.cluster_assignment]),
                priority=i,
            )
            bam_view.add_tracks(track)
    finally:
        os.chdir(curdir)

def make_home_url(pathname):
    names = {
        '/woldlab/castor/home/sau/public_html/': 'http://woldlab.caltech.edu/~sau/',
        '/woldlab/loxcyc/home/diane/proj/': 'http://woldlab.caltech.edu/~diane/',
    }
    for name in names:
        if pathname.startswith(name):
            return pathname.replace(name, names[name])

def hex_to_ucsc_color(c):
    colors = [str(int(c[1:3], 16)), str(int(c[3:5], 16)), str(int(c[5:7], 16))]
    return ",".join(colors)

if __name__ == '__main__':
    main()
