import numpy as np
import time

from ..helpers import errorgenerator
from ..helpers import tools
from .. import datastructures
from .. import settings
from ..helpers import generators
from .. import properties

import logging
logger = logging.getLogger(__name__)



class IdealPCR():

    def __init__(self, settings: settings.PCRSettings = settings.PCRSettings(), **kwargs):
        self.settings = settings
        self.settings.update(**kwargs)


    def __repr__(self):
        return f"{type(self).__name__}(primers={'+'.join(self.settings.primers)}, ncycles={self.settings.n_cycles})"


    #
    # Public functions
    #

    def process(self, pool: datastructures.SeqPool):
        """ Main entry point. Performs the PCR by calling the amplify function for every cycle. """

        # set up
        t = time.time()
        logger.info(f"Starting PCR with {self.settings.n_cycles} cycles.")
        pool.simplify()

        # do sampling by volume if set
        if self.settings.template_volume:
            pool = pool.sample_by_volume(self.settings.template_volume, remove_sampled_oligos=True)

        amplifiable_pool, nonamplifiable_pool = self._find_amplifiable_sequences(pool)
        logger.info(f"Amplifiable: {amplifiable_pool.n_sequences} seq / {amplifiable_pool.n_oligos} oligos")
        logger.info(f"Non-Amplifiable: {nonamplifiable_pool.n_sequences} seq / {nonamplifiable_pool.n_oligos} oligos")

        if amplifiable_pool.n_sequences > 0:
            # simulate the per-cycle counts for all amplifiable sequences
            amplified_pool = self._amplify(amplifiable_pool)

            # add the primer sequences to both ends
            amplified_pool = generators.attach_primers_to_pool(amplified_pool, datastructures.Seq(self.settings.primers[0]), datastructures.Seq(self.settings.primers[1]))
            
            # add non-amplified oligos
            amplified_pool.combine_with(nonamplifiable_pool)
            amplified_pool.volume = self.settings.volume
            amplified_pool.is_doublestranded = True
        else:
            logger.warning("Pool cannot be amplified.")
            amplified_pool = nonamplifiable_pool

        logger.info(f"PCR finished in {time.time()-t:.2f} seconds.")
        return amplified_pool


    # 
    # Pool filtering and processing
    # 

    def _amplify(self, pool: datastructures.SeqPool):

        # get initial copy counts
        final_copies = np.array(list(pool.counts()), dtype=np.int64)

        # get list of sequences and efficiencies
        seq_list = list(pool.sequences())
        efficiencies = self._get_efficiency(seq_list)
        logger.info(f"Total of {len(seq_list)} sequences, with efficiency {100*np.mean(efficiencies):.2f}% (sd {100*np.std(efficiencies):.2f}%, min {100*np.min(efficiencies):.2f}%, max {100*np.max(efficiencies):.2f}%).")

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
        nonamplifiable_pool = datastructures.SeqPool(is_doublestranded=seqpool.is_doublestranded)

        primer1_fw = datastructures.Seq(self.settings.primers[0])
        primer2_rc = datastructures.Seq(self.settings.primers[1])
        primer1_rc = primer1_fw.reverse_complement()
        primer2_fw = primer2_rc.reverse_complement()

        for sequence, count in seqpool:

            # we first assume the sequence to be in the normal orientation and look for the first primer
            pos_fwd = sequence.find_primer(primer1_fw, min_overlap=self.settings.primer_minoverlap)

            # if the first primer is found, find the second primer in the normal orientation as well
            if pos_fwd:
                pos_rvs = sequence.rfind_primer(primer2_fw, min_overlap=self.settings.primer_minoverlap)

                # if we found the second primer as well, we are done
                if pos_rvs:
                    properties.AmplificationEfficiency.duplicate_efficiencies(sequence, [sequence[pos_fwd:pos_rvs]])
                    amplifiable_pool.add_sequence(sequence[pos_fwd:pos_rvs], count)
                    continue

            # if the forward orientation was not successful, try the reverse orientation and look for the second primer
            pos_fwd = sequence.find_primer(primer2_rc, min_overlap=self.settings.primer_minoverlap)

            # if the second primer is found, find the first primer in the reverse orientation as well
            if pos_fwd:
                pos_rvs = sequence.rfind_primer(primer1_rc, min_overlap=self.settings.primer_minoverlap)

                # if we found the first primer as well, we are done but need to reverse-complement the sequence into normal orientation
                if pos_rvs:
                    seq_rc = sequence[pos_fwd:pos_rvs].reverse_complement()
                    properties.AmplificationEfficiency.duplicate_efficiencies(sequence, [seq_rc])
                    amplifiable_pool.add_sequence(seq_rc, count)
                    continue

            # if we didn't find any primers, the sequence is not amplifiable
            nonamplifiable_pool.add_sequence(sequence, count)

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
      
    #
    # Overwrite ideal functions
    #
 
    def _get_copies(self, sequences, final_copies):
        """ Returns a SeqPool containing the copied and imperfect sequences. """

        # as first condition, use all sequences within 0.1% of the most common sequence
        threshold = float(np.max(final_copies)) / 1000
        ix_abovethreshold = np.nonzero(final_copies > threshold)[0]
        ix_belowthreshold = np.nonzero(final_copies <= threshold)[0]
        logger.debug(f"Initial condition has sequences above threshold: {len(ix_abovethreshold)}, below threshold: {len(ix_belowthreshold)}. Coverage: {100*np.sum(final_copies[ix_abovethreshold], dtype=float)/np.sum(final_copies, dtype=float):.2f}%")

        # as second condition, make sure this covers at least the specified ratio of all oligos
        if np.sum(final_copies[ix_abovethreshold])/np.sum(final_copies) < 0.95:
            logger.debug("Changing threshold to cover sufficient oligos.")
            # sort by abundance and find the cutoff sequence after which threshold is met
            ix2sort = np.argsort(final_copies)[::-1]
            ix_delimiter = np.argmax(np.cumsum(final_copies[ix2sort]) > 0.95*np.sum(final_copies)) + 1
            # update the indexes based on the cutoff
            ix_abovethreshold = ix2sort[0:ix_delimiter]
            ix_belowthreshold = ix2sort[ix_delimiter:]

        # log stats
        logger.debug(f"Sequences above threshold: {len(ix_abovethreshold)}, below threshold: {len(ix_belowthreshold)}.")
        logger.debug(f"Oligos above threshold: {np.sum(final_copies[ix_abovethreshold], dtype=float)}/{np.sum(final_copies, dtype=float)} ({100*np.sum(final_copies[ix_abovethreshold], dtype=float)/np.sum(final_copies, dtype=float):.2f}%)")

        # add all sequences which be mutated
        amplified_pool = datastructures.SeqPool()
        amplified_pool.add_sequences(map(sequences.__getitem__, ix_abovethreshold), final_copies[ix_abovethreshold])

        # mutate the sequences
        substituter = errorgenerator.Substitutions(
            mean_rate=self.settings.n_cycles*self.settings.polymerase_basesubstitutionrate/self.settings.polymerase_fidelity, 
            bias=self.settings.polymerase_substitutionbias,
            length_bias=None,
            read_bias=None,
            repeat_bias=None,
            error_generation_error_coverage=self.settings.error_generation_error_coverage,
            error_generation_max_coverage_threshold=self.settings.error_generation_max_coverage_threshold,
        )
        amplified_pool = substituter.mutate_pool(amplified_pool)

        # add all sequences which are below the threshold
        amplified_pool.add_sequences(map(sequences.__getitem__, ix_belowthreshold), final_copies[ix_belowthreshold])
        return amplified_pool


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