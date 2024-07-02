import configargparse
import pathlib
import subprocess
import time

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

import dt4dds
import dt4dds.analysis.errormapper
import dt4dds.analysis.config
dt4dds.default_logging()

# 
# individual sub-routines 
# 

def errors(folderpath, options):
    """ Analyze the errors in the reads, using a pre-generated mapping file. """
    start = time.time()

    # get the read file
    bam_file = folderpath / dt4dds.analysis.config.bam_file

    # if using force, regenerate the mapping file
    if options.force:
        logger.warning('Forcing re-generation of mapping file.')
        mapping(folderpath, options, only_coverage=False)

    # if it doesn't exist, generate the mapping file
    elif not bam_file.exists():
        logger.warning('Generating the mapping file, as it is not present.')
        mapping(folderpath, options, only_coverage=False)

    # if it still doesn't exist, raise an error
    if not bam_file.exists():
        logger.critical(f'Mapping file not found at {bam_file}.')
        raise FileNotFoundError(f'Mapping file not found at {bam_file}.')

    # run the error analysis
    dt4dds.analysis.errormapper.errors(
        bam_file,
        reference_file=pathlib.Path(folderpath) / dt4dds.analysis.config.design_file
    )
    logger.info(f'Finished error analysis after {int(time.time()-start)}s.')


def mapping(folderpath, options, only_coverage=False):
    """ Perform the mapping of the input files to a sam file using bbmap. """
    start = time.time()
    if options.paired:
        logger.warning('Using paired reads for mapping.')
        files = f'in1={folderpath / dt4dds.analysis.config.input_file1} in2={folderpath / dt4dds.analysis.config.input_file2}'
        interleaved = 'int=t'
    else:
        logger.warning('Using only forward reads for mapping.')
        files = f'in1={folderpath / dt4dds.analysis.config.input_file1}'
        interleaved = 'int=f'

    if options.cut_adapters:
        logger.warning('Cutting adapters from reads. This may cause false-positive errors at the ends of reads due to over-trimming.')
        cutting = " ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo"
    else:
        cutting = ""

    if only_coverage:
        logger.warning('Generating only coverage statistics.')
        output = f"scafstats={folderpath / dt4dds.analysis.config.coverage_file}"
    else:
        output = f"mappedonly out={folderpath / dt4dds.analysis.config.bam_file} scafstats={folderpath / dt4dds.analysis.config.output_dir / dt4dds.analysis.config.coverage_file}"

    if options.map_local:
        logger.warning('Using local mapping with bbmap.')
        aln_mode = "local=t"
    else:
        aln_mode = "local=f"

    # remove old mapping file
    try:
        (folderpath / dt4dds.analysis.config.bam_file).unlink()
        logger.warning('Removed old mapping file.')
    except FileNotFoundError:
        pass

    # build commands
    command1 = f"{dt4dds.analysis.config.bbmap_path}bbduk.sh -Xmx{4*dt4dds.analysis.config.n_workers}g {files} samplerate={dt4dds.analysis.config.align_subsampling} out=stdout.fq t={dt4dds.analysis.config.n_workers}{cutting}"
    command2 = f"{dt4dds.analysis.config.bbmap_path}bbmap.sh -Xmx{4*dt4dds.analysis.config.n_workers}g in=stdin.fq t={dt4dds.analysis.config.n_workers} ref={folderpath / dt4dds.analysis.config.design_file} {interleaved} {aln_mode} nodisk vslow k=8 maxindel=200 minratio=0.1 {output}"
    logger.info(f'Running mapping with bbmap: {command1} | {command2}')
    
    # run process
    output = subprocess.run(f'{command1} | {command2}', shell=True, capture_output=True, text=True)
    print(output.stderr)

    # write output to file
    with open(folderpath / dt4dds.analysis.config.output_dir / 'log_bbmap.out', 'w') as f:
        f.write(output.stdout)
        f.write(output.stderr)

    # check for errors
    if output.returncode != 0 or 'Exception' in output.stderr:
        logger.critical('Error running bbmap: {}'.format(output.stderr))
        raise RuntimeError('Error running bbmap: {}'.format(output.stderr))
    logger.info(f'Finished mapping with bbmap after {int(time.time()-start)}s.')


# config file parser which allows for selection of default files
def config_file_open_func(filename, mode='r'):
    """ Custom config file opener which allows for selection of default files. """
    if mode == 'w': return open(filename, 'w')
    # if we are reading a config file, first get the default files
    default_files = {p.stem: p.resolve() for p in (pathlib.Path(__file__).parent / 'analysis-configs').glob('*.yaml')}
    # try to get the filename from the defaults, if its not there, it wasn't a default file and we just open it as usual
    return open(default_files.get(filename, filename), 'r')
    

# create parser
parser = configargparse.ArgParser(
    config_file_parser_class=configargparse.YAMLConfigFileParser,
    config_file_open_func=config_file_open_func,
    description='Analysis of experiments'
)

# required arguments
parser.add('folderpath', type=str, help='Path to the experiment folder')
parser.add('-p', '--paired', help='Use both forwards and reverse reads for analysis', action='store_true')
parser.add('-f', '--force', help='Force re-generation of the alignment map from raw reads', action='store_true')
parser.add('-l', '--local', help='Perform a local alignment against the reference', action='store_true')
parser.add('-rm', '--readmetrics', help='Generate and save read metrics', action='store_true')
parser.add('-bm', '--breakmetrics', help='Generate and save breakage metrics', action='store_true')
parser.add('-ra', '--cut_adapters', help='Cut adapters from reads before alignment', action='store_true')
parser.add('-ml', '--map_local', help='Enable local mapping and alignment with bbmap', action='store_true')
parser.add('-dir', '--categorize_by_direction', help='Aggregate reads by read direction rather than read number', action='store_true')

# config file
parser.add('-c', '--config', is_config_file=True, help='Path to a config file')

# system parameters
parser.add('-w', '--n_workers', type=int, help='Number of workers/cores to use for analysis', default=dt4dds.analysis.config.n_workers)

# analysis parameters
parser.add('-s', '--align_subsampling', type=float, help='Subsampling rate for alignment', default=dt4dds.analysis.config.align_subsampling)
parser.add('-n', '--n_reads', type=int, help='Number of reads to use for analysis', default=dt4dds.analysis.config.n_reads)
parser.add('--min_match_similarity', type=float, help='Minimum similarity for a match', default=dt4dds.analysis.config.min_match_similarity)
parser.add('--medium_match_similarity', type=float, help='Medium similarity for a match', default=dt4dds.analysis.config.medium_match_similarity)
parser.add('--min_align_length', type=int, help='Minimum length of an alignment', default=dt4dds.analysis.config.min_align_length)
parser.add('--trim_ends', type=str, help='Bases to trim from the ends of the reads, as part of the sequencing adapters', nargs='*', default=dt4dds.analysis.config.trim_ends)
parser.add('--trim_adapter', type=str, help='(Sub-)Sequence to trim from the end of the reads, i.e. sequencing adapter', default=dt4dds.analysis.config.trim_adapter)


# get and update parameters
options = parser.parse_args()
for key, value in options.__dict__.items():
    if hasattr(dt4dds.analysis.config, key):
        setattr(dt4dds.analysis.config, key, value)
logger.info("Using parameters: {}".format(options))

def run(folderpath):
    # perform the requested analysis type
    try:
        (folderpath / dt4dds.analysis.config.output_dir).mkdir(exist_ok=True)
        errors(folderpath, options)
    except:
        logger.exception('Error running analysis on {}'.format(folderpath))


def generate_command(folderpath):
    # create a config file with the parameters
    mode = 'local' if options.local else 'global'
    config_file_path = folderpath / dt4dds.analysis.config.output_dir / f"{mode}_{dt4dds.analysis.config.params_file}"
    (folderpath / dt4dds.analysis.config.output_dir).mkdir(exist_ok=True)
    parser.write_config_file(options, [str(config_file_path),])
    
    # compile the command and submit it to the cluster
    command = f'dt4dds-analysis {folderpath} -c {config_file_path}'
    return command


# single entry point
def single():
    # get the experiment folder
    folderpath = pathlib.Path(options.folderpath).resolve()
    if not folderpath.exists():
        raise FileNotFoundError('Experiment folder not found: {}'.format(folderpath))
    
    # run the analysis
    run(folderpath)

    
# batch entry point
def batch():
    # get the parent folder
    folderpath = pathlib.Path(options.folderpath).resolve()
    if not folderpath.exists():
        raise FileNotFoundError('Experiment folder not found: {}'.format(folderpath))
    
    # run the analysis for each subfolder
    for subfolder in folderpath.glob('*/'):
        logger.info('Running analysis on {}'.format(subfolder))
        run(subfolder)


def cluster_single():
    # get the experiment folder
    folderpath = pathlib.Path(options.folderpath).resolve()
    if not folderpath.exists():
        raise FileNotFoundError('Experiment folder not found: {}'.format(folderpath))
    
    # run the analysis
    dt4dds_command = generate_command(folderpath)
    logger.info(f'Generated analysis command: {dt4dds_command}')
    command = dt4dds.analysis.config.cluster_submission_command(folderpath, options, dt4dds_command)
    logger.info(f'Submitting analysis command to cluster: {command}')
    subprocess.run(command, shell=True, capture_output=True, text=True)


def cluster_batch():
    # get the parent folder
    folderpath = pathlib.Path(options.folderpath).resolve()
    if not folderpath.exists():
        raise FileNotFoundError('Experiment folder not found: {}'.format(folderpath))
    
    # run the analysis for each subfolder
    for subfolder in folderpath.glob('*/'):
        if not subfolder.is_dir(): continue
        dt4dds_command = generate_command(subfolder)
        logger.info(f'Generated analysis command: {dt4dds_command}')
        command = dt4dds.analysis.config.cluster_submission_command(subfolder, options, dt4dds_command)
        logger.info(f'Submitting analysis command to cluster: {command}')
        subprocess.run(command, shell=True, capture_output=True, text=True)
        time.sleep(0.2)