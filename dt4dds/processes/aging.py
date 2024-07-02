import numpy as np
import tqdm.auto
import time

from ..helpers import errorgenerator
from ..helpers import tools
from .. import datastructures
from .. import settings
from ..helpers import config

import logging
logger = logging.getLogger(__name__)


def _breakage_model(sequence, count, cut_rate, base_probabilities={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}):

    # holds the new oligos
    new_oligos = []

    # generate number of breaks
    init_len = len(sequence)
    n_breaks = tools.rng.binomial(init_len, cut_rate, size=count)

    # get position breakage probability for this sequence
    base_probabilities['_'] = 0.0
    probs = np.array([base_probabilities[base] for base in sequence if base in base_probabilities])
    probs = probs / np.sum(probs)

    # generate break positions
    pos_list = iter(tools.rng.choice(init_len, p=probs, size=sum(n_breaks)))

    # for each oligo, generate the broken oligos based on positions of the breaks
    for n_breaks_i in n_breaks:

        # short-circuit if no breaks
        if n_breaks_i == 0:
            new_oligos.append([sequence,])
            continue

        # get the break positions and sort them
        break_positions = np.fromiter(pos_list, count=n_breaks_i, dtype='H')
        break_positions.sort()

        # will hold all broken oligos
        oligos = []

        # add the first broken oligo
        if break_positions[0] != 0:
            oligo = sequence[:break_positions[0]] + '_'
            oligos.append(oligo)

        # add the middle broken oligos
        for i in range(1, n_breaks_i):
            start = break_positions[i-1] + 1 # +1 to skip the previous break
            end = break_positions[i]
            oligo = sequence[start:end] + '_'
            oligos.append(oligo)

        # add the last broken oligo
        if break_positions[-1] != init_len-1:
            oligo = sequence[break_positions[-1]+1:]
            oligos.append(oligo)

        # add the oligos to the list
        new_oligos.append(oligos)

    return new_oligos




class IdealAging():


    def __init__(self, settings: settings.AgingSettings = settings.AgingSettings(), n_halflives=None, intact_ratio=None, decayed_ratio=None, **kwargs):

        self.settings = settings
        self.settings.update(**kwargs)

        # calculate the decay rates
        if n_halflives is not None:
            logger.info("Using n_halflives to calculate decay progress.")
            self.n_halflives = n_halflives
        elif intact_ratio is not None:
            logger.info("Using intact_ratio to calculate decay progress.")
            self.n_halflives = -np.log2(intact_ratio)
        elif decayed_ratio is not None:
            logger.info("Using decayed_ratio to calculate decay progress.")
            self.n_halflives = -np.log2(1-decayed_ratio)
        else:
            raise ValueError("Must specify either n_halflives, intact_ratio, or decayed_ratio.")


    def __repr__(self):
        return f"{type(self).__name__}(n_halflives={self.n_halflives:.2f}, intact_ratio={self.intact_ratio:.4f}, decayed_ratio={self.decayed_ratio:.4f})"


    @property
    def intact_ratio(self):
        return 2**(-self.n_halflives)
    
    @property
    def decayed_ratio(self):
        return 1 - self.intact_ratio


    #
    # Public functions
    #

    def process(self, pool: datastructures.SeqPool):
        """ Main entry point. Performs decay. """

        # set up
        t = time.time()
        logger.info(f"Starting decay with {self}.")
        pool.simplify()

        # introduce errors and decay
        error_pool = self._errors(pool)
        decayed_pool = self._decay(error_pool)

        logger.info(f"Aging finished in {time.time()-t:.2f}s.")
        return decayed_pool


    # 
    # Decay-related
    # 

    def _decay(self, pool):
        if not self.settings.breakage_enabled:
            return pool.sample_by_factor(self.intact_ratio, remove_sampled_oligos=False)
        else:
            return self._breakage(pool)
        

    def _errors(self, pool):
        return pool


    # 
    # Breakage-related
    # 

    def _breakage(self, pool):
        """ Breaks oligos into pieces. """

        # calculate cut rate
        if self.settings.breakage_rate_override:
            cut_rate = self.settings.breakage_rate_override
        else:
            cut_rate = 1 - (self.intact_ratio)**(1/self.settings.breakage_ref_length)
        logger.info(f"Breakage cut rate: {cut_rate:.4f} per nt")
        
        # break oligos to form a new pool
        broken_pool = datastructures.SeqPool(is_doublestranded=False)
        for sequence, count in tqdm.tqdm(pool, total=pool.n_sequences, desc="Breaking", unit="sequences", disable=not config.show_progressbars, leave=False):

            # add the original sequence
            self._breakage_build_pool(broken_pool, sequence, count, cut_rate)

            # add the reverse complement if the pool is double-stranded
            if pool.is_doublestranded:
                self._breakage_build_pool(broken_pool, sequence.reverse_complement(), count, cut_rate)

        return broken_pool


    def _breakage_build_pool(self, broken_pool, sequence, count, cut_rate):
    
        if count > self.settings.breakage_coverage_threshold:
            count_factor = count // self.settings.breakage_coverage_threshold
            overflow = count % self.settings.breakage_coverage_threshold
            count = self.settings.breakage_coverage_threshold
        else:
            count_factor = 1
            overflow = 0

        # break oligos
        new_oligos = _breakage_model(str(sequence), count, cut_rate, base_probabilities=self.settings.breakage_base_probabilities)
        for broken_oligos in new_oligos:
            broken_pool.add_single_sequences(broken_oligos, count=count_factor)
        
        # add overflow
        if overflow > 0:
            ixs = tools.rng.integers(len(new_oligos), size=overflow)
            for ix in ixs:
                broken_pool.add_single_sequences(new_oligos[ix], count=1)





class Aging(IdealAging):

    #
    # Overwrite ideal functions
    #

    def _errors(self, sample_pool):
        """  """
        # initiate the error generators with their settings
        error_generators = []
        error_generators.append(errorgenerator.Substitutions(
            mean_rate=self.settings.substitution_rate*self.n_halflives, 
            bias=self.settings.substitution_bias,
            read_bias=self.settings.substitution_read_bias,
            repeat_bias=self.settings.substitution_repeat_bias,
            length_bias=self.settings.substitution_length_bias,
            error_generation_error_coverage=self.settings.error_generation_error_coverage,
            error_generation_max_coverage_threshold=self.settings.error_generation_max_coverage_threshold,
        ))
        error_generators.append(errorgenerator.Deletions(
            mean_rate=self.settings.deletion_rate*self.n_halflives, 
            bias=self.settings.deletion_bias, 
            read_bias=self.settings.deletion_read_bias,
            repeat_bias=self.settings.deletion_repeat_bias,
            length_bias=self.settings.deletion_length_bias,
            error_generation_error_coverage=self.settings.error_generation_error_coverage,
            error_generation_max_coverage_threshold=self.settings.error_generation_max_coverage_threshold,
        ))
        error_generators.append(errorgenerator.Insertions(
            mean_rate=self.settings.insertion_rate*self.n_halflives, 
            bias=self.settings.insertion_bias,
            read_bias=self.settings.insertion_read_bias,
            repeat_bias=self.settings.insertion_repeat_bias,
            length_bias=self.settings.insertion_length_bias,
            error_generation_error_coverage=self.settings.error_generation_error_coverage,
            error_generation_max_coverage_threshold=self.settings.error_generation_max_coverage_threshold,
        ))

        # apply the error generators
        for error_gen in error_generators:
            sample_pool = error_gen.mutate_pool(sample_pool)

        return sample_pool
