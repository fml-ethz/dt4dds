import random
import tqdm
import time
import itertools

from ..helpers import errorgenerator
from .. import datastructures
from .. import settings
from ..helpers import tools
from ..helpers import config

import logging
logger = logging.getLogger(__name__)




class IdealSBSSequencing():

    def __init__(self, settings: settings.SequencingSettings = settings.SequencingSettings(), **kwargs):
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
    # Public functions
    #

    def process(self, pool: datastructures.SeqPool):
        """ Main entry point. Performs the Sequencing. """

        # start timer
        t = time.time()
        logger.info(f"Starting sequencing with outputs to {self.settings.output_directory}.")

        # select reads based on primers and total number of reads
        subpool = self._select_sequences(pool)

        # separate and mutate the reads
        reads = self._generate_reads(subpool)

        # save the reads
        tools.seqlist_to_fastq(reads[0], f"{self.settings.output_directory}/R1.fq.gz", use_gzip=True)
        if self.paired_end: 
            tools.seqlist_to_fastq(reads[1], f"{self.settings.output_directory}/R2.fq.gz", use_gzip=True)

        logger.info(f"Sequencing completed in {time.time()-t:.2f} seconds.")

    # 
    # Sequence filtering
    # 

    def _select_sequences(self, seqpool: datastructures.SeqPool):
        """ Limit to sequences that can be sequenced and select up to the specified limit. """

        # only sequences with adequate primers can be bridge-amplified and sequenced
        sequenceable_pool = self._find_sequenceable_sequences(seqpool)

        # sample oligos to get expected number of reads
        if sequenceable_pool.n_oligos > self.settings.n_reads:
            sequenceable_pool = sequenceable_pool.sample_by_counts(self.settings.n_reads, remove_sampled_oligos=False, accurate=True)
            logger.info(f"Selected {sequenceable_pool.n_oligos} oligos of {sequenceable_pool.n_sequences} sequences for sequencing.")
        elif sequenceable_pool.n_oligos == 0:
            logger.exception("Unable to sequence, no sequence-able oligos found.")
            raise ValueError("No sequence-able oligos found.")
        else:
            logger.warning(f"Only {sequenceable_pool.n_oligos} oligos available for total of {self.settings.n_reads} reads. Continuing.")

        return sequenceable_pool


    def _find_sequenceable_sequences(self, seqpool: datastructures.SeqPool):
        """ Find all sequences that can be sequenced, e.g. have the corresponding primer sequences, and order them so that forwards read is in the normal orientation. """

        sequencable_pool = datastructures.SeqPool()

        pool_iterator = tqdm.tqdm(
            seqpool, 
            desc="Selecting", 
            total=seqpool.n_sequences,
            disable=not config.show_progressbars,
        )

        primer1_fw = datastructures.Seq(self.settings.primers_read[0])
        primer2_rc = datastructures.Seq(self.settings.primers_read[1])
        primer1_rc = primer1_fw.reverse_complement()
        primer2_fw = primer2_rc.reverse_complement()

        
        for sequence, count in pool_iterator:
            
            # we first assume the sequence to be in the normal orientation and look for the first primer
            pos_fwd = sequence.find_primer(primer1_fw, min_overlap=self.settings.primer_minoverlap)

            # if the first primer is found, find the second primer in the normal orientation as well
            if pos_fwd:
                pos_rvs = sequence.rfind_primer(primer2_fw, min_overlap=self.settings.primer_minoverlap)

                # if we found the second primer as well, we are done
                if pos_rvs:
                    sequencable_pool.add_sequence(sequence, count)
                    continue

            # if the forward orientation was not successful, try the reverse orientation and look for the second primer
            pos_fwd = sequence.find_primer(primer2_rc, min_overlap=self.settings.primer_minoverlap)

            # if the second primer is found, find the first primer in the reverse orientation as well
            if pos_fwd:
                pos_rvs = sequence.rfind_primer(primer1_rc, min_overlap=self.settings.primer_minoverlap)

                # if we found the first primer as well, we are done but need to reverse-complement the sequence into normal orientation
                if pos_rvs:
                    sequencable_pool.add_sequence(sequence.reverse_complement(), count)
                    continue

        return sequencable_pool


    # 
    # Read processing
    # 

    def _generate_reads(self, seqpool: datastructures.SeqPool):
        """ Generate reads from the pool of sequences. """
        reads_iterator = self._categorize_reads(seqpool)
        mutated_reads_iterator = self._mutate_reads(reads_iterator)
        padtrimmed_reads_iterator = self._pad_trim_reads(mutated_reads_iterator)

        reads = tqdm.tqdm(
            padtrimmed_reads_iterator, 
            desc="Sequencing", 
            total=self.settings.n_reads,
            disable=not config.show_progressbars,
        )
        # go from iterator of (fw, rv) tuples to a list of shuffled fw and rv reads
        readlist = list(reads)
        random.shuffle(readlist)
        return list(zip(*readlist))
        


    def _categorize_reads(self, pool: datastructures.SeqPool):
        ''' Separates the raw reads into forward and reverse reads. '''
        for seq, count in pool:
            read_fw = self._get_read_from_sequence(seq, self.settings.primers_read[0])
            # if we are doing paired-end sequencing, also get the reverse read
            if self.paired_end:
                seq_rc = seq.reverse_complement()
                read_rv = self._get_read_from_sequence(seq_rc, self.settings.primers_read[1])
                yield read_fw, read_rv, count
            else:
                yield read_fw, count


    def _get_read_from_sequence(self, sequence: datastructures.Seq, primer: str):
            ''' Assumes the sequences have been ordered such that the primer is at the start. '''
            # start by scanning the forward strand for the sequencing primer
            seq_start = sequence.find_primer(primer, min_overlap=self.settings.primer_minoverlap)
            if seq_start:
                return sequence.sequence[seq_start:seq_start+self.settings.read_length]

            # if we did not find anything, something went wrong
            logger.exception(f"Unable to find primer {primer} in sequence {sequence}.")
            return None


    def _mutate_reads(self, reads_iterator):
        ''' Returns the perfect reads from the pool. No mutations introduced here. '''
        # just yield from the reads_iterator, but expand each read by the count
        for read_data in reads_iterator:
            if self.paired_end:
                read_fw, read_rv, count = read_data
                for _ in range(count):
                    yield (read_fw, read_rv)
            else:
                read_fw, count = read_data
                for _ in range(count):
                    yield (read_fw,)


    def _pad_trim_reads(self, reads_iterator):
        """ Pad or trim all reads in list to full read length. """
        # create a cycle of random bases to pad with
        pad_bases = itertools.cycle(random.choices(tools.BASES, k=500))

        # define a function to pad or trim a sequence to the read length
        def padtrimmer(sequence):
            seq_l = len(sequence)
            if seq_l < self.settings.read_length:
                sequence = sequence + ''.join([next(pad_bases) for _ in range(self.settings.read_length-seq_l)])
            elif seq_l > self.settings.read_length:
                sequence = sequence[0:self.settings.read_length]
            return sequence

        # pad or trim all reads
        for read_data in reads_iterator:
            read_fw = padtrimmer(read_data[0])
            if self.paired_end:
                read_rv = padtrimmer(read_data[1])
                yield (read_fw, read_rv)
            else:
                yield (read_fw,)




class SBSSequencing(IdealSBSSequencing):

    #
    # Overwrite ideal functions
    #

    def _mutate_reads(self, reads_iterator):
        ''' Returns the mutated reads based on the specified error generators. '''
        yield from self._error_generator(self._pool_generator(reads_iterator))


    def _pool_generator(self, reads_iterator):
        """ Returns a generator that encapsulates each read into its own SeqPool for error generation. """

        for read_data in reads_iterator:
            if self.paired_end:
                read_fw, read_rv, count = read_data
                pool_fw = datastructures.SeqPool(pool_data={read_fw: count})
                pool_rv = datastructures.SeqPool(pool_data={read_rv: count})
                yield (pool_fw, pool_rv)
            else:
                read_fw, count = read_data
                pool_fw = datastructures.SeqPool(pool_data={read_fw: count})
                yield (pool_fw,)


    def _error_generator(self, pools_iterator):

        # initiate the error generators with their settings
        error_generators = []
        error_generators.append([
            errorgenerator.Substitutions(
                mean_rate=self.settings.substitution_rates[i], 
                bias=self.settings.substitution_bias[i],
                read_bias=self.settings.substitution_read_bias[i],
                repeat_bias=self.settings.substitution_repeat_bias[i],
                length_bias=self.settings.substitution_length_bias[i],
                error_generation_error_coverage=self.settings.error_generation_error_coverage,
                error_generation_max_coverage_threshold=self.settings.error_generation_max_coverage_threshold,
            ) for i in range(2 if self.paired_end else 1)])
        
        # go from iterator of (fw, rv) tuples to a tuple of iterators of fw and rv
        iters = itertools.tee(pools_iterator, len(self.read_modes))
        mutated_pools = [(pools[0] for pools in iters[0])]
        if self.paired_end: mutated_pools.append((pools[1] for pools in iters[1]))

        # sequentially apply each error generator to the corresponding pool
        for error_gen in error_generators:
            for i in range(2 if self.paired_end else 1):
                mutated_pools[i] = error_gen[i].mutate_iterator(mutated_pools[i], disable_progressbar=True)

        # yield the mutated reads, back in the original format of (fw, rv) tuples
        if self.paired_end:
            for pool_fw, pool_rv in zip(*mutated_pools):
                for seq_fw, seq_rv in zip(pool_fw.to_list(), pool_rv.to_list()):
                    yield seq_fw, seq_rv
        else:
            for pool_fw in mutated_pools[0]:
                for seq_fw in pool_fw.to_list():
                    yield (seq_fw,)