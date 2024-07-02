from dataclasses import dataclass
import random
import pysam
import functools
import Bio.SeqIO
import multiprocessing
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

from . import mapping
from . import config

def read_generator(bam_file):
    logger.info(f"Using input file {bam_file} for read generation.")
    with pysam.AlignmentFile(bam_file, "rb", threads=config.n_workers) as bamfile:
        for read in bamfile.fetch(until_eof=True):
            if read.is_mapped:
                ref_id = bamfile.get_reference_name(read.reference_id)
                yield read, ref_id


def reservoir_sampler(bam_file, n):
    logger.info(f"Reading input file {bam_file} to determine total number of reads.")
    total_reads = int(pysam.view('-c', str(bam_file)).strip("\n"))
    selected_reads = set(random.sample(range(total_reads), min(n, total_reads)))
    logger.info(f"Randomly selecting {len(selected_reads)} reads from {total_reads} total reads.")

    max_seq = max(selected_reads)
    for i, record in enumerate(read_generator(bam_file)):
        if not i % (total_reads//5):
            logger.info(f"Read {i}/{total_reads} ({100*i/total_reads:.1f}%) reads.")
        if i in selected_reads:
            yield record
        if i == max_seq:
            break


def initial_sampler(bam_file, n):
    logger.warning(f"Selecting only the first {n} reads from all reads.")
    for record in read_generator(bam_file):
        if n:
            n -= 1
            yield record
        else:
            break






def read_reference(input_file):
    d = dict()
    with open(input_file) as f:
        for record in Bio.SeqIO.parse(f, 'fasta'):
            d[record.id] = record

    logger.info(f"Generated {len(d)} reference sequences from {input_file}.")
    return d


def reverse_complement_references(reference_dict):
    d = dict()
    for id, ref in reference_dict.items():
        d[id] = ref[:]
        d[id].seq = ref.seq.reverse_complement()
    logger.info(f"Reverse complemented the {len(d)} reference sequences.")
    return d



def ensure_no_duplicate_references(reference_dict):
    result_dict = {}
    known_seqs = set()
    for ref_id, ref in reference_dict.items():
        if ref.seq not in known_seqs:
            known_seqs.add(ref.seq)
            result_dict[ref_id] = ref

    if len(result_dict) == len(reference_dict):
        logger.info("No duplicates in reference sequence list.")
    else:
        logger.warning(f"Removed {len(reference_dict)-len(result_dict)} duplicates from reference sequence list.")
    return result_dict







def _sum_stats(mappings_list):
    stats = [mapping.ErrorMap(), mapping.ErrorMap()]

    if config.categorize_by_direction:
        categorize = lambda x: x.read_direction-1
    else:
        categorize = lambda x: x.read_number-1

    for map in mappings_list:
        stats[categorize(map)] += mapping.ErrorMap.from_mapping(map)
    
    return stats





def _save_error_stats(mappings_list, output_files):

    if not mappings_list:
        logger.info(f"No data to write to error stats {mappings_list}.")
        return
    
    n_per_worker = -(-len(mappings_list) // config.n_workers)
    separated_mappings = [mappings_list[i:i+n_per_worker] for i in range(0, len(mappings_list), n_per_worker)]

    with multiprocessing.Pool(config.n_workers) as pool:
        results = pool.map(_sum_stats, separated_mappings)

    stats = results[0]
    for res in results[1:]:
        stats[0] += res[0]
        stats[1] += res[1]

    for i, file in enumerate(output_files):
        stats[i].save(file)

    logger.info(f"Wrote error statistics to files {', '.join([str(x) for x in output_files])}.")



def _save_error_maps(mappings_list, output_files, limit=1000):

    if not mappings_list:
        logger.info(f"No data to write to error maps {', '.join([str(x) for x in output_files])}.")
        return
    
    
    if config.categorize_by_direction:
        categorize = lambda x: x.read_direction-1
    else:
        categorize = lambda x: x.read_number-1

    fhs = [open(file, "w") for file in output_files]
    written = [0, 0]
    for map in mappings_list:
        if written[categorize(map)] < limit:
            fhs[categorize(map)].write(map.to_comparison_string())
            written[categorize(map)] += 1
    
    logger.info(f"Wrote error maps to files {', '.join([str(x) for x in output_files])}, with a limit of {limit}.")




def save_results(folderpath, mappings_dict):

    mapping_files_args = []
    stat_files_args = []
    mode = 'local' if config.local else 'global'

    for name, mappings in mappings_dict.items():
        mapping_files_args.append((
            mappings, 
            [folderpath / f"{d}.{mode}.{name}.{config.errormap_suffix}" for d in ('fw', 'rv')],
        ))
        stat_files_args.append((
            mappings, 
            [folderpath / f"{d}.{mode}.{name}.{config.stats_suffix}" for d in ('fw', 'rv')],
        ))

    mapping_files_fun = functools.partial(_save_error_maps, **{'limit': config.mapping_limit})
    stats_files_fun = functools.partial(_save_error_stats)

    for map_args, stat_args in zip(mapping_files_args, stat_files_args):
        mapping_files_fun(*map_args)
        stats_files_fun(*stat_args)



def save_metrics(folderpath, mappings_dict, function_name, suffix):
    mode = 'local' if config.local else 'global'
    with open(folderpath / f'{mode}.{suffix}', "w") as f: 
        # extract the necessary columns from the keys of the first element
        columns = list(getattr(next(iter(mappings_dict.values()))[0], function_name).keys())
        columns.append('category')
        f.write(','.join(columns) + '\n')
        
        # loop through all mappings and compile the metrics
        for name, maps in mappings_dict.items():
            for read in maps:
                metrics = list(getattr(read, function_name).values())
                metrics.append(name)
                f.write(f"{','.join([str(x) for x in metrics])}\n")

    logger.info(f"Wrote {function_name} to file {folderpath / f'{mode}.{suffix}'}.")