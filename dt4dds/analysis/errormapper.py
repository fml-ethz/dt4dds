import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

from . import alignment
from . import fileio
from . import config



def errors(bam_file, reference_file):

    logger.info(f"Running error analysis on {bam_file} with reference {reference_file}.")

    # read references and reverse-complement them
    references_fw = fileio.ensure_no_duplicate_references(fileio.read_reference(reference_file))
    references_rv = fileio.reverse_complement_references(references_fw)

    # build a generator for the reads
    reads = fileio.reservoir_sampler(bam_file, config.n_reads)

    # map the reads and categorize them based on similarity
    mappings = alignment.align_reads(reads, references_fw, references_rv)
    maps_categorized = alignment.categorize_mappings(mappings)

    # save the results
    if config.categorize_by_direction: logger.warning("Categorizing by direction instead of read number.")
    fileio.save_results(bam_file.parent / config.output_dir, maps_categorized)

    # save read metrics if required
    if config.readmetrics:
        logger.info(f"Saving read metrics for {bam_file}.")
        fileio.save_metrics(bam_file.parent / config.output_dir, maps_categorized, 'readmetrics', config.readmetrics_suffix)

    # save break metrics if required
    if config.breakmetrics:
        logger.info(f"Saving break metrics for {bam_file}.")
        fileio.save_metrics(bam_file.parent / config.output_dir, maps_categorized, 'breakmetrics', config.breakmetrics_suffix)