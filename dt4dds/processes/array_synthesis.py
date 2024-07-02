import numpy as np
import random
import itertools

from .. import datastructures
from .. import settings
from .. import properties
from ..helpers import tools
from ..helpers import errorgenerator

import logging
logger = logging.getLogger(__name__)


class IdealArraySynthesis():

    def __init__(self, settings: settings.SynthesisSettings = settings.SynthesisSettings(), **kwargs):

        self.settings = settings
        self.settings.update(**kwargs)

        self.pool = datastructures.SeqPool(is_doublestranded=False)
    

    def __repr__(self):
        return f"{type(self).__name__}()"


    #
    # Public functions
    #

    def process(self, sequence_list: list):
        """ Main entry point. Performs the Sequencing. """

        # convert mean oligo scale in fmol to molecules
        mean_oligo_count = round(self.settings.per_oligo_scale * 1e-15 * tools.NA)

        # create efficiencies for the sequences
        properties.AmplificationEfficiency.get_efficiencies(sequence_list)

        # create the base pool
        self._synthesize(sequence_list, mean_oligo_count)


    # 
    # Pool manipulation
    # 

    def _synthesize(self, sequence_list: list, mean_oligo_count: int):
        """"""
        self.pool.add_sequences(sequence_list, [mean_oligo_count]*len(sequence_list))


    def sample_by_factor(self, dil_factor, volume=None):
        """"""
        # determine if scale-down is required, if so perform it
        coverage = dil_factor*self.pool.n_oligos/self.pool.n_sequences
        if (scale_factor := int(coverage/1000)) > 1:
            dil_factor /= scale_factor

        logger.info(f"Coverage: {coverage:0.1f}, scale_factor: {scale_factor}, dil_factor: {dil_factor}.")

        sample_pool = self.pool.sample_by_factor(
            dil_factor, 
            remove_sampled_oligos=True
        )
        sample_pool = self._process_sampled_oligos(sample_pool)
        sample_pool.volume = volume
        sample_pool.is_doublestranded = False

        # if scale-down was performed, reverse it
        if scale_factor > 1:
            for seq in sample_pool.sequences():
                sample_pool[seq] *= scale_factor

        return sample_pool


    def sample_by_counts(self, counts, volume=None):
        """"""
        dil_factor = counts/self.pool.n_oligos
        return self.sample_by_factor(dil_factor, volume=volume)


    def sample_by_mass(self, mass, volume=None):
        """"""
        dil_factor = mass/self.pool.mass
        return self.sample_by_factor(dil_factor, volume=volume)


    def sample_by_concentration(self, concentration, volume):
        """"""
        mass = concentration*volume
        return self.sample_by_mass(mass, volume=volume)
   

    def _process_sampled_oligos(self, sample_pool):
        """"""
        return sample_pool






class ArraySynthesis(IdealArraySynthesis):

    #
    # Overwrite ideal functions
    #

    def _synthesize(self, sequence_list: list, mean_oligo_count: int):
        """"""
        # generate oligo count population
        oligo_counts = self._get_count_distribution(mean_oligo_count, len(sequence_list))
        self.pool.add_sequences(sequence_list, oligo_counts)


    def _get_count_distribution(self, mean_oligo_count: int, n_sequences: int):
        """ """
        
        if not hasattr(tools.rng, self.settings.oligo_distribution_type):
            raise NotImplementedError(f"A distribution called {self.settings.oligo_distribution_type} is not implemented.")
        dist = getattr(tools.rng, self.settings.oligo_distribution_type)

        oligo_counts = float(mean_oligo_count)*dist(**self.settings.oligo_distribution_params, size=n_sequences)
        oligo_counts = np.rint(oligo_counts.clip(min = 0))

        return oligo_counts
    

    def _adjust_oligo_lengths(self, pool: datastructures.SeqPool):
        """  """
        cut_off_lengths = itertools.cycle(random.choices(
            np.arange(len(self.settings.length_distribution)),
            weights=np.array(self.settings.length_distribution)/sum(self.settings.length_distribution),
            k=500
        ))

        sample_pool = datastructures.SeqPool(is_doublestranded=False)
        for seq, count in pool:
            # get substraction lengths
            new_length = np.fromiter(cut_off_lengths, count=count, dtype=np.int64)
            ix_nonfull_length = np.where(new_length)[0]
            n_full_lengths = count - len(ix_nonfull_length)

            # add full length sequences and short-cut if all are full length
            if n_full_lengths: sample_pool.add_sequence(seq, n_full_lengths)
            if n_full_lengths == count: continue

            # fill in the oligos that were cut
            oligos = []

            # generate start and end positions for the remaining sequences
            starts, ends = self._generate_length_positions(new_length[ix_nonfull_length])
            for start, end in zip(starts, ends):
                try:
                    if end > 0:
                        oligos.append(seq[start:-end])
                    else:
                        oligos.append(seq[start:])
                except Exception as e:
                    logger.exception(f"Error in sequence {seq} with start {start} and end {end}.")
                    raise e

            # copy the efficiency property for all new oligos
            properties.AmplificationEfficiency.duplicate_efficiencies(seq, oligos)

            # add the oligos to the pool
            sample_pool.add_single_sequences(oligos)

        return sample_pool


    def _generate_length_positions(self, lengths: np.ndarray):

        if self.settings.length_distribution_location == 'start':
            starts = np.zeros_like(lengths)
            ends = lengths

        elif self.settings.length_distribution_location == 'end':
            starts = lengths
            ends = np.zeros_like(lengths)

        elif self.settings.length_distribution_location == 'even':
            starts = tools.rng.integers(0, np.max(lengths)+1, size=len(lengths))
            starts = np.mod(starts, np.array(lengths)+1)
            ends = lengths - starts

        return starts, ends



    def _process_sampled_oligos(self, sample_pool: datastructures.SeqPool):
        """  """
        # initiate the error generators with their settings
        error_generators = []
        error_generators.append(errorgenerator.Substitutions(
            mean_rate=self.settings.substitution_rate, 
            bias=self.settings.substitution_bias,
            read_bias=self.settings.substitution_read_bias,
            repeat_bias=self.settings.substitution_repeat_bias,
            length_bias=self.settings.substitution_length_bias,
            error_generation_error_coverage=self.settings.error_generation_error_coverage,
            error_generation_max_coverage_threshold=self.settings.error_generation_max_coverage_threshold,
        ))
        error_generators.append(errorgenerator.Deletions(
            mean_rate=self.settings.deletion_rate, 
            bias=self.settings.deletion_bias, 
            read_bias=self.settings.deletion_read_bias,
            repeat_bias=self.settings.deletion_repeat_bias,
            length_bias=self.settings.deletion_length_bias,
            error_generation_error_coverage=self.settings.error_generation_error_coverage,
            error_generation_max_coverage_threshold=self.settings.error_generation_max_coverage_threshold,
        ))
        error_generators.append(errorgenerator.Insertions(
            mean_rate=self.settings.insertion_rate, 
            bias=self.settings.insertion_bias,
            read_bias=self.settings.insertion_read_bias,
            repeat_bias=self.settings.insertion_repeat_bias,
            length_bias=self.settings.insertion_length_bias,
            error_generation_error_coverage=self.settings.error_generation_error_coverage,
            error_generation_max_coverage_threshold=self.settings.error_generation_max_coverage_threshold,
        ))

        # reverse all sequences, because synthesis is in 3' -> 5' direction
        reverse_pool = datastructures.SeqPool(is_doublestranded=False)
        for seq, count in sample_pool:
            seq_r = seq.reverse()
            # copy the efficiency property for the reverse sequence
            properties.AmplificationEfficiency.overwrite_sequence(seq, seq_r)
            reverse_pool.add_sequence(seq_r, count)

        # mutate the sequences
        for error_gen in error_generators:
            reverse_pool = error_gen.mutate_pool(reverse_pool)

        # reverse the sequences back to normal orientation
        sample_pool = datastructures.SeqPool(is_doublestranded=False)
        for seq, count in reverse_pool:
            seq_r = seq.reverse()
            # copy the efficiency property for the now new forward sequence
            properties.AmplificationEfficiency.overwrite_sequence(seq, seq_r)
            sample_pool.add_sequence(seq_r, count)
        
        # implement the oligo length distribution, if specified
        if self.settings.length_distribution is not None:
            logger.info("Adjusting oligo lengths.")
            sample_pool = self._adjust_oligo_lengths(sample_pool)

        return sample_pool