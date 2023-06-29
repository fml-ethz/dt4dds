from dataclasses import dataclass
import random
import mimetypes
import functools
import gzip
import Bio.SeqIO
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

from . import mapping


def reservoir_sampler(iterable, n):
    reservoir = []
    for t, item in enumerate(iterable):
        if t < n:
            reservoir.append(item)
        else:
            m = random.randint(0, t)
            if m < n:
                reservoir[m] = item
    logger.info(f"Reservoir sampler sampled {len(reservoir)} items.")
    return reservoir


def initial_sampler(iterable, n):
    reservoir = []
    for t, item in enumerate(iterable):
        if t < n:
            reservoir.append(item)
        else:
            break
    logger.info(f"Initial sampler sampled {len(reservoir)} items.")
    return reservoir





def single_read_generator(input_file):
    logger.info(f"Single reads are being generated from {input_file}.")
    encoding = mimetypes.guess_type(input_file)[1]  # uses file extension
    open_call = functools.partial(gzip.open, input_file, mode='rt') if encoding == 'gzip' else functools.partial(open, input_file, mode='r')

    with open_call() as f:
        for record in Bio.SeqIO.parse(f, 'fastq'):
            yield record


def paired_read_generator(input_file_fw, input_file_rv):
    logger.info(f"Paired reads are being generated from {input_file_fw} (fw) and {input_file_rv} (rv).")
    open_calls = []
    for input_file in (input_file_fw, input_file_rv):
        encoding = mimetypes.guess_type(input_file)[1] #uses file extension
        open_calls.append(
            functools.partial(gzip.open, input_file, mode='rt') if encoding == 'gzip' 
            else functools.partial(open, input_file, mode='r')
        )

    with open_calls[0]() as f_fw, open_calls[1]() as f_rv:
            for record_fw, record_rv in zip(Bio.SeqIO.parse(f_fw, 'fastq'), Bio.SeqIO.parse(f_rv, 'fastq')):
                yield (record_fw, record_rv)




def read_reference(input_file):
    d = dict()
    with open(input_file) as f:
        for record in Bio.SeqIO.parse(f, 'fasta'):
            d[record.id] = record

    logger.info(f"Generated reference sequences from {input_file}.")
    return d


def reverse_complement_references(reference_dict):
    d = dict()
    for id, ref in reference_dict.items():
        d[id] = ref[:]
        d[id].seq = ref.seq.reverse_complement()
    logger.info(f"Reverse complemented the reference sequences.")
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



def _save_error_stats(mappings, output_files, paired=False):

    if not mappings:
        logger.info(f"No data to write to error stats {', '.join(output_files)}.")
        return

    stats = [mapping.ErrorMap() for _ in output_files]

    for map in mappings:
        if paired:
            for imap in map:
                stats[imap.read_number-1] += mapping.ErrorMap.from_mapping(imap)
        else:
            stats[map.read_number-1] += mapping.ErrorMap.from_mapping(map)

    for i, file in enumerate(output_files):
        stats[i].save(file)

    logger.info(f"Wrote error statistics to files {', '.join(output_files)}.")



def _save_error_maps(mappings, output_files, limit=10000, paired=False):

    if not mappings:
        logger.info(f"No data to write to error maps {', '.join(output_files)}.")
        return

    fhs = [open(file, "w") for file in output_files]

    for i, map in enumerate(mappings):
        if paired:
            for imap in map:
                fhs[imap.read_number-1].write(imap.to_comparison_string())
        else:
            fhs[map.read_number-1].write(map.to_comparison_string())

        if limit and i == limit:
            break
    
    logger.info(f"Wrote error maps to files {', '.join(output_files)}, with a limit of {limit}.")




def save_results(mappings_dict, output_files, map_limit=10000, paired=False):

    mapping_files_args = []
    stat_files_args = []
    paired_addendum = ".paired" if paired else ""

    mappings_files_suffix = "errormap"
    for name, mappings in mappings_dict.items():
        mapping_files_args.append(
            (
                mappings, 
                [f"{file}{paired_addendum}.{name}.{mappings_files_suffix}" for file in output_files],
            )
        )

    stats_files_suffix = "stats"
    for name, mappings in mappings_dict.items():
        stat_files_args.append(
            (
                mappings, 
                [f"{file}{paired_addendum}.{name}.{stats_files_suffix}" for file in output_files],
            )
        )

    mapping_files_fun = functools.partial(_save_error_maps, **{'limit': map_limit, 'paired': paired})
    stats_files_fun = functools.partial(_save_error_stats, **{'paired': paired})

    for map_args, stat_args in zip(mapping_files_args, stat_files_args):
        mapping_files_fun(*map_args)
        stats_files_fun(*stat_args)


def _save_seqerror_matrix(mappings_dict, output_file):

    stats_by_sequence = {}
    for refid, list_of_maps in mappings_dict.items():
        if list_of_maps: 
            stats = mapping.ErrorMap()
            for map in list_of_maps:
                stats += mapping.ErrorMap.from_mapping(map)
            stats_by_sequence[refid] = stats

    # compile overall error rates
    dataframes_by_parameter = {'n_reads': [], 'n_bases': [], 'n_substitutions': [], 'n_deletions': [], 'n_insertions': []}

    for refid, stats in stats_by_sequence.items():
        for parameter in dataframes_by_parameter.keys():
            dataframes_by_parameter[parameter].append(getattr(stats, parameter))
        
    df = pd.DataFrame(dataframes_by_parameter, index=list(stats_by_sequence.keys()))
    df.to_csv(f"{output_file}.overview_by_sequence.csv", index_label="refid")


    # compile error rate by position
    dataframes_by_error = {'substitutions': [], 'deletions': [], 'insertions': []}

    # collect per-seq stats
    for refid, stats in stats_by_sequence.items():
        for errortype in dataframes_by_error.keys():
            df = pd.DataFrame.from_dict(getattr(stats, f"n_{errortype}_by_refposition"), orient="index", columns=[refid])
            df[refid] /= stats.n_reads
            dataframes_by_error[errortype].append(df)

    # compile into one dataframe
    for errortype, dfs in dataframes_by_error.items():
        all_df = pd.DataFrame().join(dfs, how='outer', sort=True)  
        all_df.fillna(0, inplace=True)
        all_df.T.to_csv(f"{output_file}.{errortype}_by_sequence.csv", index_label="refid")


def save_seqerror_matrix(mappings, output_files):
    for m, f in zip(mappings, output_files):
        _save_seqerror_matrix(m, f)
