import numpy as np
import tqdm.auto
import multiprocessing

from ..helpers import config
from .. import datastructures
from .. import settings
from ..helpers import tools
from ..helpers import errorgenerator
from ..helpers import config
from .. import properties

from ..helpers.step import Step

import logging
logger = logging.getLogger(__name__)




class IdealArraySynthesis(Step):

    #
    # Basic functions
    #

    def __init__(self, settings: settings.SynthesisSettings = settings.SynthesisSettings(), **kwargs):
        super().__init__()
        
        
        self.settings = settings
        self.settings.update(**kwargs)

        self.pool = datastructures.SeqPool(is_doublestranded=False)
    



    #
    # Internal processing functions
    #

    def process(self, sequence_list: list):
        """ Main entry point. Performs the Sequencing. """

        # set up
        self._run_pre_process_hooks()

        # convert mean oligo scale in fmol to molecules
        mean_oligo_count = round(self.settings.per_oligo_scale * 1e-15 * tools.NA)

        # get the pool
        self._synthesize(sequence_list, mean_oligo_count)

        # finish up
        self._run_post_process_hooks()




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
    # Basic functions
    #

    def __init__(self, settings: settings.SynthesisSettings  = settings.SynthesisSettings(), **kwargs):

        super().__init__(settings, **kwargs)

        # initiate the error generators with their settings
        self.substituter = errorgenerator.Substitutions(
            self.settings.substitution_rate, 
            self.settings.substitution_bias
        )
        self.deleter = errorgenerator.DeletionEvents(
            self.settings.deletion_rate, 
            self.settings.deletion_bias, 
            read_bias = self.settings.deletion_read_bias,
            read_bias_params = self.settings.deletion_read_bias_params,
            repeat_bias = self.settings.deletion_repeat_bias,
            repeat_bias_params = self.settings.deletion_repeat_bias_params,
            length_bias = self.settings.deletion_length_bias,
            length_bias_params = self.settings.deletion_length_bias_params,
        )
        self.inserter = errorgenerator.Insertions(
            self.settings.insertion_rate, 
            self.settings.insertion_bias
        )

        

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



    def _process_sampled_oligos(self, sample_pool: datastructures.SeqPool):

        sampled_pool = datastructures.SeqPool()

        if config.enable_multiprocessing:
        
            with multiprocessing.Pool(processes=config.n_processes) as processpool:
                sample_pools = processpool.imap(
                    self._process_sampled_oligos_single, 
                    zip(sample_pool.sequences(), sample_pool.counts()),
                    chunksize=min(500, sample_pool.n_sequences//config.n_processes)
                )

                pools_iterator = tqdm.auto.tqdm(
                    sample_pools, 
                    desc="Mutating", 
                    unit="sequences",
                    total=sample_pool.n_sequences,
                    disable=not config.show_progressbars,
                    leave=False,
                )

                # add new sequences to pool
                for orig_seq, pool in zip(sample_pool.sequences(), pools_iterator):
                    sampled_pool.add_sequences(pool.sequences(), pool.counts())
                    properties.AmplificationEfficiency.duplicate_efficiencies(orig_seq, pool.sequences())

        else:

            sample_pools = map(self._process_sampled_oligos_single, zip(sample_pool.sequences(), sample_pool.counts()))

            pools_iterator = tqdm.auto.tqdm(
                sample_pools, 
                desc="Mutating", 
                unit="sequences",
                total=sample_pool.n_sequences,
                disable=not config.show_progressbars,
                leave=False,
            )

            # add new sequences to pool
            for orig_seq, pool in zip(sample_pool.sequences(), pools_iterator):
                sampled_pool.add_sequences(pool.sequences(), pool.counts())
                properties.AmplificationEfficiency.duplicate_efficiencies(orig_seq, pool.sequences())

        return sampled_pool



    def _process_sampled_oligos_single(self, arguments):

        sequence, count = arguments
        sample_pool = datastructures.SeqPool(pool_data = {sequence.reverse(): count})

        for mutation in (self.substituter, self.inserter, self.deleter):
            new_pool = datastructures.SeqPool()

            for sequence, count in sample_pool:
                new_pool.combine_with(mutation(sequence, count))
            sample_pool = new_pool

        new_pool = datastructures.SeqPool()
        for sequence, count in sample_pool:
            new_pool.add_sequence(sequence.reverse(), count)
        
        return new_pool