#!/usr/bin/python3
from argparse import ArgumentParser
import os
import pandas
from generate_combined_transcript_C1 import ASOF_RUN17_experiment_files, ASOF_RUN17_library_files
from woldrnaseq import models

def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('-o', '--output', help='output directory')
    parser.add_argument('--mode', default=None, choices=[
        'customtrack',
        'trackhub',
        'merge_paper_wiggles',
        'paper_median_coverage',
        'check_bedgraphs',
        'localize_tsvs',
        'paper_as_single_experiment_tsv',
        'paper_as_cluster_experiment_tsv',
    ])
    args = parser.parse_args(cmdline)

    experiment_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_experiment_files.split() ]
    library_files = [ os.path.expanduser(x.strip()) for x in ASOF_RUN17_library_files.split() ]

    experiments = models.load_experiments(experiment_files)
    libraries = models.load_library_tables(library_files)

    to_include = read_peng_20180710_cluster_memberships()
    #print('{} cells to include'.format(len(to_include)))

    if args.mode == 'customtrack':
        make_custom_tracks()
    elif args.mode == 'trackhub':
        make_trackhub()
    elif args.mode == 'merge_paper_wiggles':
        merge_paper_wiggles(to_include, libraries)
    elif args.mode == 'paper_median_coverage':
        make_paper_median_coverage(to_include, libraries, args.output)
    elif args.mode == 'check_bedgraphs':
        check_bedgraphs(to_include, libraries)
    elif args.mode == 'localize_tsvs':
        localize_tsvs(experiments, libraries, args.output)
    elif args.mode == 'paper_as_single_experiment_tsv':
        paper920_as_single_experiment_tsv(to_include, args.output)
    elif args.mode == 'paper_as_cluster_experiment_tsv':
        paper920_as_cluster_experiment_tsv(to_include, args.output)
    else:
        parser.error('Did you want to pick an operation mode?')


def make_paper_median_coverage(to_include, libraries, output_dir):
    gtf_path = os.path.expanduser('~/proj/genome/mm10-M4-male/gencode.vM4-tRNAs-ERCC.gff')
    source_type = None
    gene_type_filter = None

    triplet = 'mm10-M4-male'
    suffix = '-' + triplet + '_uniq.bw'
    count = 0
    coverages = []

    lines = ['universe=vanilla',
             'output=logs/coverage.$(process).out',
             'error=logs/coverage.$(process).out',
             'log=logs/coverage.$(process).log',
             'executable=/usr/bin/python3',
             'request_memory=4G',
             '']

    for cell_id, wigfile, cluster in find_paper_wigfiles(to_include, libraries, suffix):
        cell_id = cell_id.replace('_mm10', '').replace('_clean', '')
        target_dir = os.path.abspath(os.path.join(output_dir, cell_id))
        if not os.path.exists(target_dir):
            os.mkdir(target_dir)
        target_name = cell_id + '_' + triplet + '.coverage'
        target_pathname = os.path.join(target_dir, target_name)
        command = 'arguments="/woldlab/loxcyc/home/diane/proj/GeorgiScripts/gene_coverage_wig_gtf.py --gtf {gtf} -o {target} --gene-type protein_coding --print-list {wigname}"'.format(
            gtf=gtf_path,
            target=target_pathname,
            wigname=wigfile,
        )
        lines.append(command)
        lines.append('queue')

    lines.append('')
    condor_file = os.path.join(output_dir, 'coverage.condor')
    with open(condor_file, 'wt') as outstream:
        outstream.write('\n'.join(lines))


def merge_paper_wiggles(to_include, libraries):
    for suffix in ["-mm10-M4-male_all.bw", "-mm10-M4-male_uniq.bw"]:
        clustered_wigs = {}
        for cell_id, wigfile, cluster in find_paper_wigfiles(to_include, libraries, suffix):
            clustered_wigs.setdefault(cluster, []).append(wigfile)

        for cluster in clustered_wigs:
            bigwig_name = 'C1_peng_20180710_cluster_bigwigs/' + cluster.replace(' ', '_') + suffix
            args = ['python3', 'merge_bw.py', '-o', os.path.abspath(bigwig_name)]
            args.extend(clustered_wigs[cluster])
            print(' '.join(args))


def find_paper_wigfiles(to_include, libraries, suffix):
    """
    """
    found = 0
    missing = 0
    for cell_metadata in to_include:
        cell_id = cell_metadata['cell_id']
        cluster = cell_metadata['cluster_name']

        analysis_dir = libraries.loc[cell_id]['analysis_dir']

        wigdir = get_wigdir(analysis_dir)
        wigname = os.path.join(wigdir, cell_id + suffix)
        wigfound = os.path.exists(wigname)
        if wigfound:
            found += 1
            yield cell_id, wigname, cluster
        else:
            missing += 1
            print(cell_id, cluster, wigname)


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


def check_bedgraphs(to_include, libraries):
    not_found = 0
    for library_id, row in libraries.iterrows():
        filename = 'Signal.Unique.str1.out.bg'
        #filename = library_id + '-mm10-M4-male_genome.bam'
        pathname = os.path.join(row.analysis_dir, filename)
        if not os.path.exists(pathname):
            not_found += 1
            print('Missing: {}'.format(pathname))
    print('missing {}'.format(not_found))


def sanitize_library_suffix(x):
    return x.replace('_mm10', '').replace('_clean', '')

def localize_tsvs(experiments, libraries, target_dir):

    experiment_filename = os.path.join(target_dir, 'sequencing-experiment.tsv')
    if not os.path.exists(experiment_filename):
        clean_replicates = []
        for reps in experiments['replicates']:
            clean_replicates.append(','.join([sanitize_library_suffix(x) for x in reps]))

        experiments['replicates'] = clean_replicates
        experiments[['replicates']].to_csv(experiment_filename, sep='\t')
    else:
        print(experiment_filename, 'exists')
    libraries_filename = os.path.join(target_dir, 'library.tsv')
    cols = ['analysis_dir', 'genome', 'annotation', 'sex', 'read_1', 'reference_prefix']
    libraries.index = [sanitize(x) for x in libraries.index]
    libraries['analysis_dir'] = [sanitize_library_suffix(os.path.split(x)[1]) for x in libraries['analysis_dir']]
    if not os.path.exists(libraries_filename):
        libraries.index.name = 'library_id'
        libraries[cols].to_csv(libraries_filename, sep='\t')
    else:
        print(libraries_filename, 'already exists')


def paper920_as_single_experiment_tsv(to_include, target_dir):
    experiment_filename = os.path.join(target_dir, 'one-experiment.tsv')
    if not os.path.exists(experiment_filename):
        clean_replicates = []
        for metadata in to_include:
            clean_replicates.append(sanitize_library_suffix(metadata['cell_id']))

        experiments = pandas.DataFrame({
            'experiments': ['920_cells'],
            'replicates': ','.join(clean_replicates)
        })
        experiments.set_index('experiments', inplace=True)
        print('total replicates', len(clean_replicates))
        print(experiments.shape)
        experiments.to_csv(experiment_filename, sep='\t')
    else:
        print(experiment_filename, 'exists')


def paper920_as_cluster_experiment_tsv(to_include, target_dir):
    experiment_filename = os.path.join(target_dir, 'cluster-experiment.tsv')
    if not os.path.exists(experiment_filename):
        clean_replicates = {}
        for metadata in to_include:
            cluster = metadata['cluster_name'].replace(' ', '_')
            cell_id = sanitize_library_suffix(metadata['cell_id'])
            clean_replicates.setdefault(cluster, []).append(cell_id)

        for rep in clean_replicates:
            clean_replicates[rep] = ','.join(clean_replicates[rep])

        experiments = pandas.DataFrame({
            'experiments': list(clean_replicates.keys()),
            'replicates': [clean_replicates[x] for x in clean_replicates]
        })
        experiments.set_index('experiments', inplace=True)
        print('total replicates', len(clean_replicates))
        print(experiments.shape)
        experiments.to_csv(experiment_filename, sep='\t')
    else:
        print(experiment_filename, 'exists')


def read_peng_20180710_cluster_memberships():
    cluster_membership = []
    colors = load_violin_plot_colors()
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
            cluster_name = ' '.join(rest_values[3:])
            color = tuple([float(x) for x in rest_values[:3]])
            rgb = tuple(colors[colors['cell type'] == cluster_name]['rgb float'].values[0])
            cluster_membership.append({
                'cell_id': cell_id,
                'value': value,
                'color': rgb,
                'cluster_name': cluster_name,
            })
        return cluster_membership

def load_violin_plot_colors():
    def parse_rgb_float(rgb):
        return [float(x) for x in rgb[1:-1].split(', ')]
    def parse_rgb_decimal(rgb):
        return [int(x) for x in rgb[1:-1].split(', ')]
    data = pandas.read_csv(
        'violin-paper/peng-violin-plot-colors.csv',
        converters={
            'rgb decimal': parse_rgb_decimal,
            'rgb float': parse_rgb_float,
        })
    return data

def make_custom_tracks():
    data = read_peng_20180710_cluster_memberships()
    df = pandas.DataFrame(data, columns=['cell_id', 'cluster_name', 'value', 'color'])
    df['rgb'] = df['color'].apply(lambda channel: ','.join([str(int(x * 255)) for x in channel]))

    colors = df[['cluster_name', 'value', 'color', 'rgb']].drop_duplicates()


    template = 'track type=bigWig name={name} description={cluster_name} visibility=full color={rgb} bigDataUrl={url}'
    for suffix in ['-mm10-M4-male_all.bw', '-mm10-M4-male_uniq.bw']:
        for i, row in colors.iterrows():
            name = row.value + suffix
            print(template.format(
                name=name.replace(' ', '_'),
                cluster_name=row.cluster_name.replace(' ', '_'),
                rgb=row.rgb,
                url='http://woldlab.caltech.edu/~diane/C1_peng_20180710_cluster_bigwigs/' + name))

def make_trackhub():
    import trackhub
    data = read_peng_20180710_cluster_memberships()
    df = pandas.DataFrame(data, columns=['cell_id', 'cluster_name', 'value', 'color'])
    df['rgb'] = df['color'].apply(lambda channel: ','.join([str(int(x * 255)) for x in channel]))

    colors = df[['cluster_name', 'value', 'color', 'rgb']].drop_duplicates()

    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name="C1_peng_20180710_cluster",
        short_label="C1_peng_20180710_cluster",
        genome="mm10",
        email="diane@caltech.edu")

    subgroups = [
        trackhub.SubGroupDefinition(
            name='multiread',
            label='multiread',
            mapping={
                'all': 'all_reads',
                'uniq': 'unique_only',
            })]

    composite = trackhub.CompositeTrack(
        name='composite',
        short_label='bigwigs',
        dimensions='dimX=multiread',
        sortOrder='multiread',
        visibility='full',
        tracktype='bigWig',
    )
    composite.add_subgroups(subgroups)
    trackdb.add_tracks(composite)

    signal_view = trackhub.ViewTrack(
        name='signal',
        view='signal',
        visibility='full',
        tracktype='bigWig')
    composite.add_view(signal_view)

    subgroup_map = {
        '-mm10-M4-male_all.bw': 'all',
        '-mm10-M4-male_uniq.bw': 'uniq'
    }
    for suffix in ['-mm10-M4-male_all.bw', '-mm10-M4-male_uniq.bw']:
        for i, row in colors.iterrows():
            name = row.value + subgroup_map[suffix]
            url = 'http://woldlab.caltech.edu/~diane/C1_peng_20180710_cluster_bigwigs/' + row.value + suffix
            #url = '../' + row.value + suffix
            track = trackhub.Track(
                name=trackhub.helpers.sanitize(name),
                long_label=row.cluster_name + ' ' + subgroup_map[suffix],
                url=url,
                #source=os.path.join('C1_peng_20180710_cluster_bigwigs/', name),
                visibility='full',
                tracktype='bigWig',
                subgroups={'multiread': subgroup_map[suffix]},
                color=row.rgb)
            signal_view.add_tracks(track)

    print(trackdb)
    trackhub.upload.upload_hub(
        hub=hub,
        host='localhost',
        remote_dir='/woldlab/loxcyc/home/diane/public_html/C1_peng_20180710_cluster_bigwigs/')


if __name__ == '__main__':
    main()
