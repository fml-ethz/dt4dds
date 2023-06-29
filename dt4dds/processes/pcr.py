import gc
import numpy as np
import tqdm.auto
import multiprocessing
import itertools

from ..helpers import errorgenerator
from ..helpers import tools
from .. import datastructures
from .. import settings
from ..helpers import config
from ..helpers import generators
from .. import properties

from ..helpers.step import Step

import logging
logger = logging.getLogger(__name__)



class IdealPCR(Step):

    #
    # Basic functions
    #

    def __init__(self, settings: settings.PCRSettings = settings.PCRSettings(), **kwargs):
        super().__init__()
        
        self.settings = settings
        self.settings.update(**kwargs)

        self.primer_fwd = datastructures.Seq(self.settings.primers[0])
        self.primer_rvs = datastructures.Seq(self.settings.primers[1])



    def __repr__(self):
        return f"{type(self).__name__}(primers={'+'.join(self.settings.primers)}, ncycles={self.settings.n_cycles})"



    #
    # Internal processing functions
    #

    def process(self, pool: datastructures.SeqPool):
        """ Main entry point. Performs the PCR by calling the amplify function for every cycle. """

        # set up
        pool.simplify()
        self._run_pre_process_hooks(pool)

        # do sampling by volume if set
        if self.settings.template_volume:
            pool = pool.sample_by_volume(self.settings.template_volume, remove_sampled_oligos=True)

        amplifiable_pool, nonamplifiable_pool = self._find_amplifiable_sequences(pool)
        logger.info(f"Amplifiable: {amplifiable_pool.n_sequences} seq / {amplifiable_pool.n_oligos} oligos")
        logger.info(f"Non-Amplifiable: {nonamplifiable_pool.n_sequences} seq / {nonamplifiable_pool.n_oligos} oligos")

        gc.collect()

        if amplifiable_pool.n_sequences > 0:
            # simulate the per-cycle counts for all amplifiable sequences
            amplified_pool = self._amplify(amplifiable_pool)

            # add the primer sequences to both ends
            amplified_pool = generators.attach_primers_to_pool(amplified_pool, self.primer_fwd, self.primer_rvs)

            # add non-amplified oligos
            amplified_pool.combine_with(nonamplifiable_pool)
            amplified_pool.volume = self.settings.volume
            amplified_pool.is_doublestranded = True

        else:
            logger.warning("Pool cannot be amplified.")
            amplified_pool = nonamplifiable_pool

        # finish up
        self._run_post_process_hooks(amplified_pool)  
        return amplified_pool



    def _amplify(self, pool: datastructures.SeqPool):

        # get initial copy counts
        final_copies = np.array(list(pool.counts()), dtype=np.int64)

        # get list of sequences and efficiencies
        seq_list = list(pool.sequences())
        efficiencies = self._get_efficiency(seq_list)
        logger.info(f"Total of {len(seq_list)} sequences, with efficiency {100*np.mean(efficiencies):.2f}% (sd {100*np.std(efficiencies):.2f}%).")

        # go through each cycle and amplify
        for cycle in range(self.settings.n_cycles):

            # do not create copy on first amplification of single-stranded DNA
            if not pool.is_doublestranded:
                pool.is_doublestranded = True
                logger.info("Pool is single-stranded, skipping first amplification.")

                # pool stays at the same counts, no amplification
                continue
 
            # generate new counts for the sequences
            final_copies += self._get_copycount(final_copies, efficiencies)

        amplified_pool = self._get_copies(
            seq_list,  
            final_copies
        )
        return amplified_pool




    def _find_amplifiable_sequences(self, seqpool: datastructures.SeqPool):

        amplifiable_pool = datastructures.SeqPool(is_doublestranded=seqpool.is_doublestranded)
        nonamplifiable_pool = datastructures.SeqPool()

        for sequence, count in seqpool:

            # identify position of forwards and backwards primers
            pos_fwd = sequence.find_primer(self.primer_fwd, min_overlap=self.settings.primer_minoverlap)
            
            # if there are no primers there's no amplification, otherwise the sequence is amplified
            if not pos_fwd: 
                nonamplifiable_pool.add_sequence(sequence, count)
                continue

            # we look at the reverse complement
            sequence_revcomp = sequence.reverse_complement()
            pos_rvs = sequence_revcomp.find_primer(self.primer_rvs, min_overlap=self.settings.primer_minoverlap)
            if not pos_rvs: 
                nonamplifiable_pool.add_sequence(sequence, count)
                continue

            # from here on the sequence has all primer requirements
            amp_seq = sequence[pos_fwd:-pos_rvs]
            amplifiable_pool.add_sequence(amp_seq, count)

        return amplifiable_pool, nonamplifiable_pool


    #
    # PCR-relevant functions
    #


    def _get_copies(self, sequences, final_copies):
        """ Returns a SeqPool containing the copied, perfect sequences. """
        copypool = datastructures.SeqPool()
        copypool.add_sequences(sequences, final_copies)
        return copypool



    def _get_copycount(self, counts, efficiencies):
        """  """
        return np.rint(counts*efficiencies)



    def _get_efficiency(self, sequences):
        """  """
        return np.full(len(sequences), self.settings.efficiency_mean)




class PCR(IdealPCR):


    def __init__(self, settings: settings.PCRSettings = settings.PCRSettings(), **kwargs):

        super().__init__(settings, **kwargs)

        
    #
    # Overwrite ideal functions
    #
 
    def _get_copies(self, sequences, final_copies):
        """ Returns a SeqPool containing the copied and imperfect sequences. """

        # as first condition, use all sequences within 0.1% of the most common sequence
        threshold = float(np.max(final_copies)) / 1000
        ix_abovethreshold = np.nonzero(final_copies > threshold)[0]
        ix_belowthreshold = np.nonzero(final_copies <= threshold)[0]
        logger.info(f"Initial condition has sequences above threshold: {len(ix_abovethreshold)}, below threshold: {len(ix_belowthreshold)}. Coverage: {100*np.sum(final_copies[ix_abovethreshold], dtype=float)/np.sum(final_copies, dtype=float):.2f}%")

        # as second condition, make sure this covers at least the specified ratio of all oligos
        if np.sum(final_copies[ix_abovethreshold])/np.sum(final_copies) < 0.95:
            logger.info("Changing threshold to cover sufficient oligos.")
            # sort by abundance and find the cutoff sequence after which threshold is met
            ix2sort = np.argsort(final_copies)[::-1]
            ix_delimiter = np.argmax(np.cumsum(final_copies[ix2sort]) > 0.95*np.sum(final_copies)) + 1
            
            # update the indexes based on the cutoff
            ix_abovethreshold = ix2sort[0:ix_delimiter]
            ix_belowthreshold = ix2sort[ix_delimiter:]

        # log stats
        logger.info(f"Sequences above threshold: {len(ix_abovethreshold)}, below threshold: {len(ix_belowthreshold)}.")
        logger.info(f"Oligos above threshold: {np.sum(final_copies[ix_abovethreshold], dtype=float)}/{np.sum(final_copies, dtype=float)} ({100*np.sum(final_copies[ix_abovethreshold], dtype=float)/np.sum(final_copies, dtype=float):.2f}%)")

        # add all sequences which will not be mutated
        amplified_pool = datastructures.SeqPool(is_doublestranded=True)
        for ix in ix_belowthreshold:
            amplified_pool.add_sequence(sequences[ix], final_copies[ix])
        
        # create views into the arrays which only yield sequences to be mutated
        seq_views = map(sequences.__getitem__, ix_abovethreshold)
        sampled_view = final_copies[ix_abovethreshold]

        # mutate the sequences
        if config.enable_multiprocessing:

            with multiprocessing.Pool(processes=config.n_processes) as pool:
                pools = pool.imap(
                    self._get_copies_single, 
                    zip(seq_views, sampled_view),
                    chunksize=min(50000, len(sequences)//config.n_processes)
                )

                pools_iterator = tqdm.auto.tqdm(
                    pools, 
                    desc="Substitution", 
                    unit="sequences",
                    total=len(sequences),
                    disable=not config.show_progressbars,
                    leave=False,
                )

                # add new sequences to pool
                for orig_seq, pool in zip(map(sequences.__getitem__, ix_abovethreshold), pools_iterator):
                    amplified_pool.add_sequences(pool.sequences(), pool.counts())
                    properties.AmplificationEfficiency.duplicate_efficiencies(orig_seq, pool.sequences())

        else:

            pools = map(self._get_copies_single, zip(seq_views, sampled_view))

            pools_iterator = tqdm.auto.tqdm(
                pools, 
                desc="Substitution", 
                unit="sequences",
                total=len(sequences),
                disable=not config.show_progressbars,
                leave=False,
            )

            # add new sequences to pool
            for orig_seq, pool in zip(map(sequences.__getitem__, ix_abovethreshold), pools_iterator):
                amplified_pool.add_sequences(pool.sequences(), pool.counts())
                properties.AmplificationEfficiency.duplicate_efficiencies(orig_seq, pool.sequences())

        return amplified_pool

        



    def _get_copies_single(self, arguments):

        sequence, copies = arguments

        substituter = errorgenerator.Substitutions(
            self.settings.n_cycles*self.settings.polymerase_basesubstitutionrate/self.settings.polymerase_fidelity, 
            self.settings.polymerase_substitutionbias,
            length_bias=None
        )
        return substituter(sequence, copies)





    def _get_copycount(self, counts, efficiencies):
        """ Returns the copy number for the provided sequences. """
        return tools.rng.binomial(counts, efficiencies)
            


    def _get_efficiency(self, sequences):
        """  """
        # get relative efficiencies for known oligos
        efficiencies = properties.AmplificationEfficiency.get_efficiencies(sequences)

        # generate absolute efficiencies
        efficiencies *= (1 + self.settings.efficiency_mean)
        efficiencies = np.clip(efficiencies-1, 0.0, 1.0)

        return efficiencies