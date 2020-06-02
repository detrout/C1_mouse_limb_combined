#!/usr/bin/python3
"""Generate trackhub for specified set...
"""
import argparse
import logging
import os
import pandas
import trackhub

from woldrnaseq import models
import woldrnaseq.make_tracks
from woldrnaseq.make_tracks import make_bigwig_track_name, make_bam_track_name

from generate_combined_transcript_C1 import ASOF_RUN17_library_files
from make_rsem_subset import load_cells_set


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
    'DarkRed': 'DarkRed',
    'deep_red': 'deep_red',
    'red': 'red',
    'black': 'black',
    'green': 'green',
    'yellow': 'yellow',
}
# convert human friendly name to machine friendly
human_cluster_mapping = {
    'DarkRed': 'darkred',
    'deep red': 'deep_red',
    'red': 'red',
    'orange': 'orange',
    'black': 'black',
    'green': 'green',
    'yellow': 'yellow',
}
colors = {
    'DarkRed': '#843c39',
    'deep red': '#843c39',
    'red': '#d62728',
    'orange': '#ffa500',
    'black': '#000000',
    'green': '#2ca02c',
    'yellow': '#fde724',
}


def main(cmdline=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--dry-run', action='store_true', default=False)
    parser.add_argument('--force', action='store_true', default=False,
                        help='ignore cached merged table')
    parser.add_argument('--hub', help="hub name", required=True)
    parser.add_argument('--short-name', help="short name", required=True)
    parser.add_argument('--long-name', help="long name", required=True)
    parser.add_argument('--signal', action="store_true", default=False,
                        help="include signal (bigwig) tracks")
    parser.add_argument('--reads', action="store_true", default=False,
                        help="include reads (bam) tracks")
    parser.add_argument(
        '--track-dir', required=True,
        help='directory relative to public_html to put track in')

    parser.add_argument('filename', nargs=1)

    args = parser.parse_args(cmdline)
    filename = args.filename[0]

    collection = asof_run17_bigwig_paths.split()
    public_dir = os.path.expanduser('~/public_html/' + args.track_dir)

    base, ext = os.path.splitext(filename)
    cache_name = base + '.csv'
    if args.force or not os.path.exists(cache_name):
        roots = [os.path.expanduser(x.strip()) for x in collection]
        df = load_cells_set(filename)
        libraries = load_asof_run17_libraries()
        libraries = pandas.merge(df, libraries, how='inner', left_index=True, right_index=True)

        libraries = add_bigwig_paths(libraries, roots)
        libraries = add_bam_paths(libraries)
        libraries.to_csv(cache_name)
    else:
        libraries = pandas.read_csv(cache_name)

    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name=args.hub,
        short_label=args.short_name,
        long_label=args.long_name,
        genome='mm10',
        email='diane@caltech.edu')

    tracks_added = False
    if args.signal:
        make_bigwig_trackhub(libraries, trackdb)
        tracks_added = True
    if args.reads:
        make_bam_trackhub(libraries, trackdb)
        tracks_added = True

    if not tracks_added:
        print("Did you want to add tracks? use --signal and/or --read")

    if not args.dry_run:
        trackhub.upload.upload_hub(
            hub=hub,
            host='localhost',
            remote_dir=public_dir)
    else:
        print(trackdb)

    print('trackhub: ' + 'http://woldlab.caltech.edu/~diane/' + args.track_dir + '/' + hub.hub + '.hub.txt')


def load_asof_run17_libraries():
    library_files = [os.path.expanduser(x.strip()) for x in ASOF_RUN17_library_files.split('\n')]
    libraries = models.load_library_tables(library_files)
    name = libraries.index.name
    libraries.index = [x.replace('_mm10', '').replace('_clean', '') for x in libraries.index]
    libraries.index.name = name

    return libraries


def add_bigwig_paths(libraries, roots):
    woldrnaseq.make_tracks.logger.setLevel(logging.ERROR)
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

            if track_name is None:
                print("Couldn't find track for {}".format(library_id))
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
            priority = i * 10 + track_index
            track = trackhub.Track(
                url=make_home_url(row[track_type]),
                name="{:04d}".format(priority) + '_' + row.analysis_name + '_' + track_type,
                visibility='full',
                tracktype='bigWig',
                subgroups={'cluster': cluster_assignment, 'multi': track_type},
                color=hex_to_ucsc_color(colors[row.cluster_assignment]),
                priority=priority,
            )
            signal_view.add_tracks(track)


def make_bam_trackhub(libraries, trackdb):
    cluster = trackhub.SubGroupDefinition(
            name='cluster',
            label='cluster',
            mapping=cluster_mapping)

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
    #os.chdir(public_dir)

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
                visibility='full',
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
        '/woldlab/castor/home/sau/flowcells/C1_e10.5_mouse_limb_run2_June20_2016':
          'http://woldlab.caltech.edu/~diane/consortium-presentation/C1_e10.5_mouse_limb_run2_June20_2016',
        '/woldlab/castor/home/sau/flowcells/H5LV3BCXY':
          'http://woldlab.caltech.edu/~diane/consortium-presentation/H5LV3BCXY/',
        '/woldlab/castor/home/sau/flowcells/C1_mouse_limb_combined_Mar_2017':
          'http://woldlab.caltech.edu/~diane/consortium-presentation/C1_mouse_limb_combined_Mar_2017/',
        '/woldlab/castor/home/sau/flowcells/H7CNTBCX2/':
          'http://woldlab.caltech.edu/~diane/consortium-presentation/H7CNTBCX2/',
        '/woldlab/castor/home/sau/flowcells/HFNLTBCX2/':
          'http://woldlab.caltech.edu/~diane/consortium-presentation/HFNLTBCX2/',
        '/woldlab/castor/home/sau/flowcells/HFNLNBCX2':
          'http://woldlab.caltech.edu/~diane/consortium-presentation/HFNLNBCX2/',
        '/woldlab/castor/home/sau/flowcells/HF7NTBCX2':
          'http://woldlab.caltech.edu/~diane/consortium-presentation/HF7NTBCX2/',
        '/woldlab/castor/home/sau/flowcells/HFNYNBCX2':
          'http://woldlab.caltech.edu/~diane/consortium-presentation/HFNYNBCX2/',
    }
    for name in names:
        if pathname.startswith(name):
            return pathname.replace(name, names[name])

    raise ValueError("{} had an unrecognized prefix path".format(pathname))


def hex_to_ucsc_color(c):
    colors = [str(int(c[1:3], 16)), str(int(c[3:5], 16)), str(int(c[5:7], 16))]
    return ",".join(colors)


if __name__ == '__main__':
    main()
