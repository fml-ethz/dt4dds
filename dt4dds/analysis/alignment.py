import logging

from . import mapping
from . import config


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)



def align_reads(records, references_dict_fw, references_dict_rv):

    mode = 'local' if config.local else 'global'
    logger.warning(f"Aligning reads in mode {mode}.")
    n_reads = [[0, 0], [0, 0]]
    n_skip = [0, 0]
   
    for read, ref_id in records:
        # skip unmapped reads
        if read.is_unmapped:
            n_skip[0] += 1
            continue

        # get read direction
        if read.is_forward:
            ref = references_dict_fw[ref_id]
            read_direction = 1
        elif read.is_reverse:
            ref = references_dict_rv[ref_id]
            read_direction = 2
        else:
            n_skip[1] += 1
            continue

        # get read number
        read_number = 1
        if read.flag & 0x80: read_number = 2

        n_reads[read_direction-1][read_number-1] += 1
        yield mapping.Mapping.from_sam(read, ref, read_direction=read_direction, read_number=read_number)
    
    # print some statistics
    tot_reads = sum(n_skip) + sum([sum(n) for n in n_reads])
    logger.info(f"Aligned {tot_reads} reads, thereof {sum(n_reads[0])} ({100*sum(n_reads[0])/tot_reads:.2f}%) forward reads, {sum(n_reads[1])} ({100*sum(n_reads[1])/tot_reads:.2f}%) reverse reads, and {sum(n_skip)} ({100*sum(n_skip)/tot_reads:.2f}%) skipped reads.")
    if sum(n_skip) > 0: logger.info(f"Skipped {n_skip[0]} ({100*n_skip[0]/sum(n_skip):.2f}%) unmapped reads and {n_skip[1]} ({100*n_skip[1]/sum(n_skip):.2f}%) reads with unknown direction.")
    if sum(n_reads[0]) > 0: logger.info(f"Aligned {sum(n_reads[0])} forward reads, thereof {n_reads[0][0]} ({100*n_reads[0][0]/sum(n_reads[0]):.2f}%) from read 1, and {n_reads[0][1]} ({100*n_reads[0][1]/sum(n_reads[0]):.2f}%) from read 2.")
    if sum(n_reads[1]) > 0: logger.info(f"Aligned {sum(n_reads[1])} reverse reads, thereof {n_reads[1][0]} ({100*n_reads[1][0]/sum(n_reads[1]):.2f}%) from read 1, and {n_reads[1][1]} ({100*n_reads[1][1]/sum(n_reads[1]):.2f}%) from read 2.")




def categorize_mappings(mappings):
     
    logger.info(f"Categorizing reads, using min. match similarity of {config.min_match_similarity}, and min. medium match similarity of {config.medium_match_similarity}.")
    
    categorized_mappings = {
        'mapped_high': [],
        'mapped_low': [],
        'unmapped': [],
    }

    for mapping in mappings:

        # short-circuit if the read is below minimum length
        if len(mapping.aligned_read_sequence) < config.min_align_length:
            categorized_mappings['unmapped'].append(mapping)
            continue

        # get the similarity score of read
        score = mapping.aligned_similarity

        # this is a low scoring match which might be a random simple collision
        if score < config.min_match_similarity:
            categorized_mappings['unmapped'].append(mapping)
            continue

        # this is a medium scoring match which might be a corrupted read
        if score < config.medium_match_similarity:
            categorized_mappings['mapped_low'].append(mapping)
            continue

        # this is a good match
        categorized_mappings['mapped_high'].append(mapping)


    total_alignments = len(categorized_mappings['mapped_high']) + len(categorized_mappings['mapped_low']) + len(categorized_mappings['unmapped'])
    logger.info(f"Categorized {total_alignments} reads into {len(categorized_mappings['mapped_high'])} matches ({len(categorized_mappings['mapped_high'])/total_alignments*100:.2f}%), {len(categorized_mappings['mapped_low'])} borderline matches ({len(categorized_mappings['mapped_low'])/total_alignments*100:.2f}%), and {len(categorized_mappings['unmapped'])} unknowns ({len(categorized_mappings['unmapped'])/total_alignments*100:.2f}%).")
    return categorized_mappings