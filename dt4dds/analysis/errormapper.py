import logging
import os
import gc

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

from . import alignment
from . import fileio



default_adapter = "AGATCGGAAGAGCG"




def single(read_file, reference_file, n_reads = 1000, ref_offset = 0, read_offset = 0, read_length = None, adapter=default_adapter):

    logger.info(f"Running on {read_file} with reference {reference_file}, using {n_reads} reads, reference offset {ref_offset}, read offset {ref_offset}, read length {read_length}, adadpter {adapter}.")

    references_fw = fileio.ensure_no_duplicate_references(fileio.read_reference(reference_file))
    references_rv = fileio.reverse_complement_references(references_fw)

    reads = fileio.reservoir_sampler(fileio.single_read_generator(read_file), n_reads)

    mappings = alignment.align_single_reads_caller(
        reads, 
        references_fw, 
        references_rv, 
        read_window=60, 
        score_cutoff=0, 
        trim_seq=adapter, 
        ref_offset=ref_offset, 
        read_offset=read_offset,
        read_length=read_length
    )

    maps_categorized = alignment.categorize_mappings(mappings)

    fileio.save_results({
        'matches': maps_categorized[0],
        'medium_scores': maps_categorized[1],
        'low_scores': maps_categorized[2],
    }, (read_file,))



def paired(fw_read_file, rv_read_file, reference_file, n_reads = 1000, ref_offset = 0, read_offset = 0, read_length = None, adapter=default_adapter):

    logger.info(f"Running on {fw_read_file} and {rv_read_file} with reference {reference_file}, using {n_reads} reads, reference offset {ref_offset}, read offset {ref_offset}, read length {read_length}, adadpter {adapter}.")
    
    references_fw = fileio.ensure_no_duplicate_references(fileio.read_reference(reference_file))
    references_rv = fileio.reverse_complement_references(references_fw)

    reads = fileio.reservoir_sampler(fileio.paired_read_generator(fw_read_file, rv_read_file), n_reads)
    gc.collect()

    mappings = alignment.align_paired_reads_caller(
        reads, 
        references_fw, 
        references_rv, 
        read_window=60, 
        score_cutoff=0, 
        trim_seq=adapter,
        ref_offset=ref_offset, 
        read_offset=read_offset,
        read_length=read_length
    )

    maps_categorized = alignment.categorize_mappings(mappings)
    gc.collect()

    # 
    # OVERALL ERROR ANALYSIS
    # 
    fileio.save_results({
        'matches': maps_categorized[0],
        'medium_scores': maps_categorized[1],
        'low_scores': maps_categorized[2],
    }, (fw_read_file, rv_read_file))






def paired_cmp(fw_read_file, rv_read_file, reference_file, n_reads = 1000, ref_offset = 0, read_offset = 0, read_length = None, adapter=default_adapter):

    logger.info(f"Running on {fw_read_file} and {rv_read_file} with reference {reference_file}, using {n_reads} reads, reference offset {ref_offset}, read offset {ref_offset}, read length {read_length}, adadpter {adapter}.")

    references_fw = fileio.ensure_no_duplicate_references(fileio.read_reference(reference_file))
    references_rv = fileio.reverse_complement_references(references_fw)

    reads = fileio.reservoir_sampler(fileio.paired_read_generator(fw_read_file, rv_read_file), n_reads)

    mappings = alignment.align_paired_reads_caller(
        reads, 
        references_fw, 
        references_rv, 
        read_window=60, 
        score_cutoff=0, 
        trim_seq=adapter,
        ref_offset=ref_offset, 
        read_offset=read_offset,
        read_length=read_length
    )

    maps_categorized = alignment.categorize_mappings(mappings, paired=True)

    maps_paired = []
    for maps in maps_categorized:
        maps_paired.extend(alignment.compare_paired_mappings(maps))
    maps_categorized_paired = alignment.categorize_mappings(maps_paired, paired=True)

    fileio.save_results({
        'matches': maps_categorized_paired[0],
        'medium_scores': maps_categorized_paired[1],
        'low_scores': maps_categorized_paired[2],
    }, (fw_read_file, rv_read_file), paired=True)





def paired_seqanalysis(fw_read_file, rv_read_file, reference_file, n_reads = 1000, ref_offset = 0, read_offset = 0, read_length = None, adapter=default_adapter):

    logger.info(f"Running on {fw_read_file} and {rv_read_file} with reference {reference_file}, using {n_reads} reads, reference offset {ref_offset}, read offset {ref_offset}, read length {read_length}, adadpter {adapter}.")
    
    references_fw = fileio.ensure_no_duplicate_references(fileio.read_reference(reference_file))
    references_rv = fileio.reverse_complement_references(references_fw)

    reads = fileio.reservoir_sampler(fileio.paired_read_generator(fw_read_file, rv_read_file), n_reads)
    gc.collect()

    mappings = alignment.align_paired_reads_caller(
        reads, 
        references_fw, 
        references_rv, 
        read_window=60, 
        score_cutoff=0, 
        trim_seq=adapter,
        ref_offset=ref_offset, 
        read_offset=read_offset,
        read_length=read_length
    )

    maps_categorized = alignment.categorize_mappings(mappings)
    gc.collect()

    # 
    # OVERALL ERROR ANALYSIS
    # 
    fileio.save_results({
        'matches': maps_categorized[0],
        'medium_scores': maps_categorized[1],
        'low_scores': maps_categorized[2],
    }, (fw_read_file, rv_read_file))

    # 
    # SEQUENCE ERROR ANALYSIS
    # 
    matches_by_sequence = alignment.seqmap_paired_mappings(maps_categorized[0])
    fileio.save_seqerror_matrix(matches_by_sequence, (fw_read_file, rv_read_file))



def single_seqanalysis(fw_read_file, reference_file, n_reads = 1000, ref_offset = 0, read_offset = 0, read_length = None, adapter=default_adapter):

    logger.info(f"Running on {fw_read_file} with reference {reference_file}, using {n_reads} reads, reference offset {ref_offset}, read offset {ref_offset}, read length {read_length}, adadpter {adapter}.")
    
    references_fw = fileio.ensure_no_duplicate_references(fileio.read_reference(reference_file))
    references_rv = fileio.reverse_complement_references(references_fw)

    reads = fileio.reservoir_sampler(fileio.single_read_generator(fw_read_file), n_reads)
    gc.collect()

    mappings = alignment.align_single_reads_caller(
        reads, 
        references_fw, 
        references_rv, 
        read_window=60, 
        score_cutoff=0, 
        trim_seq=adapter,
        ref_offset=ref_offset, 
        read_offset=read_offset,
        read_length=read_length
    )

    maps_categorized = alignment.categorize_mappings(mappings)
    gc.collect()

    # 
    # OVERALL ERROR ANALYSIS
    # 
    fileio.save_results({
        'matches': maps_categorized[0],
        'medium_scores': maps_categorized[1],
        'low_scores': maps_categorized[2],
    }, (fw_read_file,))

    # 
    # SEQUENCE ERROR ANALYSIS
    # 
    matches_by_sequence = alignment.seqmap_single_mappings(maps_categorized[0])
    fileio.save_seqerror_matrix(matches_by_sequence, (fw_read_file,))





def phix(fw_read_file, rv_read_file, max_read_length=None, n_reads = 1000):

    logger.info(f"Running on {fw_read_file} and {rv_read_file}, using {n_reads} reads with {max_read_length} maximum read length.")

    reference_file = os.path.dirname(os.path.realpath(__file__)) + "/phix.fasta"

    references_fw = fileio.ensure_no_duplicate_references(fileio.read_reference(reference_file))
    references_rv = fileio.reverse_complement_references(references_fw)
    references = {'phix_fw': references_fw['phix'], 'phix_rv': references_rv['phix']}

    ref_fw = references['phix_fw']
    ref_fw.id = 'phix_fw'
    ref_rv = references['phix_rv']
    ref_rv.id = 'phix_rv'

    reads = fileio.reservoir_sampler(fileio.paired_read_generator(fw_read_file, rv_read_file), n_reads)

    maps = alignment.align_phix_reads_caller(reads, references['phix_fw'], references['phix_rv'], max_read_length=max_read_length, read_window=30)
    maps_categorized = alignment.categorize_phix_mappings(maps)

    fileio.save_results({
        'matches': maps_categorized[0],
        'medium_scores': maps_categorized[1],
        'low_scores': maps_categorized[2],
        'ambiguous': maps_categorized[3],
    }, (fw_read_file, rv_read_file))