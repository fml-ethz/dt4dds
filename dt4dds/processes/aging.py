import gc
import numpy as np
import tqdm.auto
import multiprocessing

from ..helpers import errorgenerator
from ..helpers import tools
from .. import datastructures
from .. import settings
from ..helpers import config
from .. import properties

from ..helpers.step import Step

import logging
logger = logging.getLogger(__name__)



class IdealAging(Step):

    #
    # Basic functions
    #

    def __init__(self, settings: settings.AgingSettings = settings.AgingSettings(), **kwargs):
        super().__init__()
        
        
        self.settings = settings
        self.settings.update(**kwargs)



    def __repr__(self):
        return f"{type(self).__name__}()"



    #
    # Internal processing functions
    #

    def process(self, pool: datastructures.SeqPool):
        """ Main entry point. Performs decay. """

        # set up
        pool.simplify()
        self._run_pre_process_hooks(pool)

        if self.settings.method == 'fixed':
            decayed_pool = self._fixeddecay(pool)
        elif self.settings.method == 'arrhenius':
            decayed_pool = self._arrheniusdecay(pool)
        else:
            logger.exception(f"Unknown decay method {self.settings.method}.") 

        aged_pool = self._getagedpool(decayed_pool)

        # finish up
        self._run_post_process_hooks(aged_pool)  
        return aged_pool



    def _fixeddecay(self, pool):

        assert (self.settings.fixed_decay_ratio >= 0.0) and (self.settings.fixed_decay_ratio <= 1.0)
        intact_ratio = 1.0 - self.settings.fixed_decay_ratio
        self.n_halflives = -np.log2(intact_ratio)
        return pool.sample_by_factor(intact_ratio, remove_sampled_oligos=False)


    def _arrheniusdecay(self, pool):

        # fitted parameters for polymer storage
        A = 4.96e20
        b = 0.074

        # time in years to time in seconds
        time_s = self.settings.arrhenius_t*365*24*60*60

        # determine rate constant in 1/s
        # saturation vapor pressure of water (https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185)
        psat = 10.0**(4.65-(1435/((self.settings.arrhenius_T+273)-64.8)))
        csat = psat/(8.3145*(self.settings.arrhenius_T+273))
        self.k0 = A*csat*np.exp(-self.settings.arrhenius_Ea/(8.3145*(self.settings.arrhenius_T+273)))*(self.settings.arrhenius_RH+b)

    
        counts = np.array(list(pool.counts()))
        lengths = np.array([len(s) for s in pool.sequences()])

        # for each sequence, decay depends on sequence length
        intact_ratios = np.exp(-self.k0*lengths*time_s)
        intact_counts = tools.rng.binomial(counts, intact_ratios)
        self.n_halflives = -np.log2(np.sum(intact_counts)/np.sum(counts))

        # build new pool
        decayed_pool = datastructures.SeqPool(is_doublestranded=pool.is_doublestranded)
        decayed_pool.add_sequences(pool.sequences(), intact_counts)

        return decayed_pool


    def _getagedpool(self, pool):
        return pool




class Aging(IdealAging):


    def __init__(self, settings: settings.AgingSettings = settings.AgingSettings(), **kwargs):

        super().__init__(settings, **kwargs)



    #
    # Overwrite ideal functions
    #

    def _getagedpool(self, sample_pool):
        """  """

        # initiate the error generators with their settings
        self.mutation_functions = []
        if self.settings.substitution_rate > 0:
            self.mutation_functions.append(errorgenerator.Substitutions(
                self.settings.substitution_rate*self.n_halflives, 
                self.settings.substitution_bias
            ))
        if self.settings.deletion_rate > 0:
            self.mutation_functions.append(errorgenerator.Deletions(
                self.settings.deletion_rate*self.n_halflives, 
                self.settings.deletion_bias, 
            ))
        if self.settings.insertion_rate > 0:
            self.mutation_functions.append(errorgenerator.Insertions(
                self.settings.insertion_rate*self.n_halflives, 
                self.settings.insertion_bias
            ))

            
        sampled_pool = datastructures.SeqPool()
        
        if config.enable_multiprocessing:
            
            with multiprocessing.Pool(processes=config.n_processes) as pool:
                sample_pools = pool.imap(
                    self._getagedpool_single, 
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

            sample_pools = map(self._getagedpool_single, zip(sample_pool.sequences(), sample_pool.counts()))

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



    def _getagedpool_single(self, arguments):
        """  """

        sequence, count = arguments
        sample_pool = datastructures.SeqPool(pool_data = {sequence: count})

        for mutation in self.mutation_functions:
            new_pool = datastructures.SeqPool()

            for sequence, count in sample_pool:
                new_pool.combine_with(mutation(sequence, count))
            sample_pool = new_pool
        
        return sample_pool
