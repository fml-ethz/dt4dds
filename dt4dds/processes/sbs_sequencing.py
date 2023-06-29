import multiprocessing
import tqdm.auto

from ..helpers import errorgenerator
from .. import datastructures
from .. import settings
from ..helpers import tools
from ..helpers import config

from ..helpers.step import Step

import logging
logger = logging.getLogger(__name__)




class IdealSBSSequencing(Step):

    #
    # Basic functions
    #

    def __init__(self, settings: settings.SequencingSettings = settings.SequencingSettings(), **kwargs):
        super().__init__()        
        
        self.settings = settings
        self.settings.update(**kwargs)

        # decide on sequencing read mode
        if self.settings.read_mode.lower() == 'single-end':
            self.paired_end = False
            self.read_modes = [1]
        elif self.settings.read_mode.lower() == 'paired-end':
            self.paired_end = True
            self.read_modes = [1, 2]
        else:
            logger.exception(f"Unknown sequencing read mode '{self.settings.read_mode}'.")
            raise NotImplementedError(f"Unknown sequencing read mode '{self.settings.read_mode}'.")



    def __repr__(self):
        return f"{type(self).__name__}()"



    #
    # Internal processing functions
    #

    def process(self, pool: datastructures.SeqPool):
        """ Main entry point. Performs the Sequencing. """

        # set up
        self._run_pre_process_hooks()

        # perform Sequencing
        self._sequence(pool)

        # finish up
        self._run_post_process_hooks()   



    def _sequence(self, seqpool: datastructures.SeqPool):

        cluster_seqpool = seqpool

        # sample oligos to get expected number of reads
        if cluster_seqpool.n_oligos >= self.settings.n_reads:
            cluster_seqpool = cluster_seqpool.sample_by_counts(self.settings.n_reads, remove_sampled_oligos=False)
        elif cluster_seqpool.n_oligos == 0:
            logger.exception("Unable to sequence, no sequence-able oligos found.")
        else:
            logger.warning(f"Only {cluster_seqpool.n_oligos} oligos available for total of {self.settings.n_reads} reads. Continuing.")

        generated_reads = self._process_reads(cluster_seqpool)

        for i_read, readlist in enumerate(generated_reads):
            tools.seqlist_to_fastq(readlist, f"{self.settings.output_directory}/R{i_read+1}.fq.gz", use_gzip=True)




    def _process_reads(self, seqpool: datastructures.SeqPool):

        cluster_reads = [[], []]

        # mutate the sequences
        if config.enable_multiprocessing:

            with multiprocessing.Pool(processes=config.n_processes) as pool:
                pools = pool.imap_unordered(
                    self._read, 
                    seqpool,
                    chunksize=min(10000, seqpool.n_sequences//config.n_processes)
                )

                pools_iterator = tqdm.auto.tqdm(
                    pools, 
                    desc="Mutating", 
                    unit="sequences",
                    total=seqpool.n_sequences,
                    disable=not config.show_progressbars,
                    leave=False,
                )

                for cluster in pools_iterator:
                    for i_read, cluster_read in enumerate(cluster):
                        if cluster_read: cluster_reads[i_read].extend(cluster_read)

        else:

            pools = map(self._read, seqpool)

            pools_iterator = tqdm.auto.tqdm(
                pools, 
                desc="Mutating", 
                unit="sequences",
                total=seqpool.n_sequences,
                disable=not config.show_progressbars,
                leave=False,
            )

            for cluster in pools_iterator:
                for i_read, cluster_read in enumerate(cluster):
                    if cluster_read: cluster_reads[i_read].extend(cluster_read)


        return cluster_reads





    def _read(self, arguments):

        sequence, count = arguments

        cluster_reads = [None, None]
        reads = [sequence, sequence.reverse_complement()]

        # generate a list of the forward and if enabled also backwards real read sequences
        for i_read, read in zip(self.read_modes, reads):
            seq_start = read.find_primer(self.settings.primers_read[i_read-1], min_overlap=self.settings.primer_minoverlap)
            if seq_start:
                # get_reads introduces all mutations, then pad them to read length if required
                oligo_reads = self._get_reads(read[seq_start:seq_start+self.settings.read_length], count, i_read)
            else:
                oligo_reads = [""]*count
            cluster_reads[i_read-1] = self._pad_trim_reads(oligo_reads)

        return cluster_reads




    def _get_reads(self, oligo, count, read_number):
        """ Return a list of perfect reads from a given count of an oligo. """
        return [oligo]*count



    def _pad_trim_reads(self, readlist):
        """ Pad or trim all reads in list to full read length. """

        for i, sequence in enumerate(readlist):

            seq_l = len(sequence)

            if seq_l < self.settings.read_length:
                readlist[i] = sequence + ''.join(tools.rng.choice(tools.BASES, size=self.settings.read_length-seq_l))
            elif seq_l > self.settings.read_length:
                readlist[i] = sequence[0:self.settings.read_length]
            else:
                continue

        return readlist












class SBSSequencing(IdealSBSSequencing):


    def __init__(self, settings: settings.SequencingSettings = settings.SequencingSettings(), **kwargs):

        super().__init__(settings, **kwargs)

        self.mutation_functions = []

        if sum(self.settings.substitution_rates) > 0:
            self.mutation_functions.append([
                errorgenerator.Substitutions(
                    self.settings.substitution_rates[i], 
                    self.settings.substitution_bias[i],
                    length_bias=self.settings.substitution_length_bias,
                    length_bias_params=self.settings.substitution_length_bias_params[i],
                ) 
            for i in range(2 if self.paired_end else 1)])

        if sum(self.settings.deletion_rates) > 0:
            self.mutation_functions.append([
                errorgenerator.Deletions(
                    self.settings.deletion_rates[i], 
                    self.settings.deletion_bias[i],
                    length_bias=self.settings.deletion_length_bias,
                    length_bias_params=self.settings.deletion_length_bias_params[i],
                ) 
            for i in range(2 if self.paired_end else 1)])

        if sum(self.settings.insertion_rates) > 0:
            self.mutation_functions.append([
                errorgenerator.Insertions(
                    self.settings.insertion_rates[i], 
                    self.settings.insertion_bias[i],
                    length_bias=self.settings.insertion_length_bias,
                    length_bias_params=self.settings.insertion_length_bias_params[i],
                )
            for i in range(2 if self.paired_end else 1)])


    #
    # Overwrite ideal functions
    #
 
    
    def _get_reads(self, oligo: datastructures.Seq, count: int, read_number: int):
        """ Return a list of mutated reads from a given count of an oligo. """

        sample_pool = datastructures.SeqPool()
        sample_pool.add_sequence(oligo, count)         

        for mutation in self.mutation_functions:
            new_pool = datastructures.SeqPool()
            for sequence, count in sample_pool:
                new_pool.combine_with(mutation[read_number-1](sequence, count))
            sample_pool = new_pool
        
        return sample_pool.to_list()
