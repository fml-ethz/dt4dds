import logging
logger = logging.getLogger(__name__)



n_workers = 4
""" Number of workers/cores to use for analysis """


# 
# analysis parameters
# 

paired = False
""" Whether the reads are analysed paired-end or not """

n_reads = 10000
""" Number of reads to use for analysis """

align_subsampling = 1.0
""" Subsampling factor for the alignment """

local = False
""" Use local alignment mode for mapping, otherwise use a global alignment """

trim_ends = []
""" Bases to trim from the ends of the reads, as part of the sequencing adapters """

trim_adapter = ""
""" (Sub-)Sequence to trim from the ends of the reads, e.g. the sequencing adapter """

min_match_similarity = 0.7
""" Minimum similarity to be considered a match """

medium_match_similarity = 0.85
""" Minimum similarity to be considered a good match """

min_align_length = 20
""" Minimum length of an alignment to be considered a match """

readmetrics = False
""" Whether to generate and save read metrics """

breakmetrics = False
""" Whether to generate and save break metrics """

cut_adapters = False
""" Whether to cut adapters from reads before alignment """

map_local = False
""" Whether to use local mapping and alignment with bbmap """


# 
# input/output file parameters
# 

input_file1 = 'R1.fq.gz'
input_file2 = 'R2.fq.gz'
""" Standard names for the input raw read files """

design_file = 'design_files.fasta'
""" Standard name for the design file """

output_dir = 'analysis'
""" Standard name for the output directory """

coverage_file = 'coverage.tsv'
""" Standard name for the coverage output file """

bam_file = 'mapped.bam'
""" Standard name for the generated alignment file """

errormap_suffix = 'errormap'
""" Standard suffix for the error map files """

stats_suffix = 'stats'
""" Standard suffix for the error statistics files """

readmetrics_suffix = 'readmetrics.csv'
""" Standard suffix for the read metrics files """

breakmetrics_suffix = 'breakmetrics.csv'
""" Standard suffix for the break metrics files """

params_file = 'dt4dds_analysis.conf'
""" Standard name for the parameter file """

bbmap_path = '~/.local/bin/bbmap/'
""" Path to the bbmap folder """

mapping_limit = 1000
""" Maximum number of error maps to write per category """

categorize_by_direction = False
""" Whether to categorize error maps by read direction or read number """



def to_dict():
    """ Returns a dictionary of all the parameters in this module """
    return {k: v for k, v in globals().items() if k[0] != '_'}



def cluster_submission_command(folderpath, options, command):
    n_cores = n_workers
    mode = 'local' if options.local else 'global'
    output = folderpath / output_dir / f'log_{mode}.out'
    
    time = '00-24:00'
    mem = '6G'
    
    return f'sbatch -n {n_cores} -t {time} --job-name="{folderpath.name}" --mem-per-cpu={mem} -o {output} --wrap="{command}"'