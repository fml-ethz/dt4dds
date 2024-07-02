import time
import numpy as np
import dataclasses
import itertools

from .. import datastructures
from .. import tools

import logging
logger = logging.getLogger(__name__)





class AbstractPreparation():

    def process(self, source_pool: datastructures.SeqPool):
        """ Main entry point. Performs the sample preparation. """
        
        t = time.time()
        logger.info(f"Starting sample preparation {self}.")

        # perform purification
        prepared_pool = self._prepare(source_pool)

        logger.info(f"Sample preparation {self} finished in {time.time()-t:.2f}s.")
        return prepared_pool


    def _prepare(self, source_pool: datastructures.SeqPool):
        raise NotImplementedError("You must call a subclass of preparation.")
    

    def __repr__(self):
        return f"{type(self).__name__}()"



@dataclasses.dataclass
class AdapterLigation(AbstractPreparation):

    adapter_3prime: str = ''
    adapter_5prime: str = ''
    final_volume: float = 50
    tail_bases: list = dataclasses.field(default_factory=list)
    tail_length: int = 0
    repair_3prime: float = 0.0


    def _generate_ligated_sequence(self, sequence: str):
        """ Generates a ligated sequence using the adapters, including tailing if enabled. """
        if self.tail_length > 0:
            new_seq = "".join(map(str, [self.adapter_5prime, sequence, "".join(tools.rng.choice(self.tail_bases, size=self.tail_length)), self.adapter_3prime]))
        else:
            new_seq = "".join(map(str, [self.adapter_5prime, sequence, self.adapter_3prime]))
        return new_seq


    def _prepare(self, source_pool: datastructures.SeqPool):
        """ Performs the sample preparation. """

        new_pool = datastructures.SeqPool(is_doublestranded=True, volume=self.final_volume)

        # keep track of repaired and skipped oligos
        n_ligated, n_skipped = 0, 0

        # perform adapter ligation
        for seq, count in source_pool:
            # only ligate if 3' end is not damaged or if ligation is possible for damaged 3' ends
            if not seq.damaged_3prime:
                n_ligated += count
                new_seq = self._generate_ligated_sequence(seq.sequence)
                new_pool.add_sequence(new_seq, count)
            elif self.repair_3prime > 0:
                if tools.rng.random() < self.repair_3prime:
                    n_ligated += count
                    new_seq = self._generate_ligated_sequence(seq.sequence)
                    new_pool.add_sequence(new_seq, count)
            else:
                n_skipped += count
                new_pool.add_sequence(seq, count)

        total_oligos = n_ligated + n_skipped
        logger.info(f"Ligated {n_ligated} ({100*n_ligated/total_oligos:0.2f}%) oligos, skipped {n_skipped} ({100*n_skipped/total_oligos:0.2f}%) oligos.")
        return new_pool
    


@dataclasses.dataclass
class BeadCleanup(AbstractPreparation):

    cutoff: int = 59
    slope: float = 0.015
    final_volume: float = 50


    def _get_probability(self, lengths: np.array):
        """ Returns the retention probability for the given lengths. """
        return np.clip(self.slope * (lengths - self.cutoff), 0, 1)
    

    def _prepare(self, source_pool: datastructures.SeqPool):
        """ Performs the sample preparation. """

        new_pool = datastructures.SeqPool(is_doublestranded=source_pool.is_doublestranded, volume=self.final_volume)
        random_numbers = itertools.cycle(tools.rng.random(size=500))
        probability_by_length = self._get_probability(np.arange(1, self.cutoff+(1/self.slope)+1))

        # keep track of repaired and skipped oligos
        n_kept, n_removed = 0, 0

        for seq, count in source_pool:
            seq_length = len(seq)

            # get threshold from sequence length
            if seq_length < self.cutoff:
                # if sequence is shorter than cutoff, it is always removed
                n_removed += count
                continue
            if seq_length > self.cutoff+(1/self.slope):
                # if sequence is long, it is always retained
                n_kept += count
                new_pool.add_sequence(seq, count)
                continue

            # otherwise, use a random number to determine whether the sequence is retained based on length-dependent probability
            probability = probability_by_length[seq_length]
            n_retained_oligos = np.count_nonzero(np.fromiter(random_numbers, count=count, dtype=float) < probability)
            n_kept += n_retained_oligos
            new_pool.add_sequence(seq, n_retained_oligos)

        total_oligos = n_kept + n_removed
        logger.info(f"Kept {n_kept} ({100*n_kept/total_oligos:0.2f}%) oligos, removed {n_removed} ({100*n_removed/total_oligos:0.2f}%) oligos.")
        return new_pool
    
    


@dataclasses.dataclass
class ssAdapterLigation(AbstractPreparation):
    """ Performs ssDNA-based library preparation using adaptase, bead cleanup, and adapter ligation. """

    adapter_3prime: str = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    adapter_5prime: str = 'CTCTTTCCCTACACGACGCTCTTCCGATCT'
    final_volume: float = 20
    tail_bases: list = dataclasses.field(default_factory=lambda: ['C', 'T'])
    tail_length: int = 8
    repair_3prime: float = 0.0


    def _prepare(self, source_pool: datastructures.SeqPool):
        """ Performs the sample preparation. """

        # first step: denaturation. Only relevant for double-stranded pools
        if source_pool.is_doublestranded:
            pool = datastructures.SeqPool(is_doublestranded=False)
            for seq, count in source_pool:
                pool.add_sequence(seq, count)
                pool.add_sequence(seq.reverse_complement(), count)
        else:
            pool = source_pool

        # second step: adapter ligation with 3' adapter
        adaptase = AdapterLigation(
            adapter_3prime=self.adapter_3prime,
            final_volume=self.final_volume,
            tail_bases=self.tail_bases,
            tail_length=self.tail_length,
            repair_3prime=self.repair_3prime,
        )
        pool = adaptase.process(pool)

        # third step: bead-based cleanup
        cleanup = BeadCleanup(
            final_volume=self.final_volume,
        )
        pool = cleanup.process(pool)

        # fourth step: adapter ligation with 5' adapter
        ligation = AdapterLigation(
            adapter_5prime=self.adapter_5prime,
            final_volume=self.final_volume,
            tail_length=0,
            repair_3prime=self.repair_3prime,
        )
        pool = ligation.process(pool)
        return pool
    


@dataclasses.dataclass
class dsAdapterLigation(AbstractPreparation):
    """ Performs dsDNA-based library preparation using adapter ligation and bead cleanup. """

    adapter_3prime: str = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    adapter_5prime: str = 'CTCTTTCCCTACACGACGCTCTTCCGATCT'
    final_volume: float = 20
    tail_bases: list = dataclasses.field(default_factory=lambda: [])
    tail_length: int = 0
    repair_3prime: float = 0.0


    def _prepare(self, source_pool: datastructures.SeqPool):
        """ Performs the sample preparation. """

        # we first filter based on intact 3' ends and add the reverse complement of each sequence
        pool = datastructures.SeqPool(is_doublestranded=False, volume=self.final_volume)
        for seq, count in source_pool:
            if not seq.damaged_3prime:
                pool.add_sequence(seq.sequence, count)
                pool.add_sequence(seq.reverse_complement().sequence, count)
            elif self.repair_3prime > 0:
                if tools.rng.random() < self.repair_3prime:
                    pool.add_sequence(seq.sequence, count)
                    pool.add_sequence(seq.reverse_complement().sequence, count)

        # first step: adapter ligation with 3' adapter
        ligation = AdapterLigation(
            adapter_3prime=self.adapter_3prime,
            adapter_5prime=self.adapter_5prime,
            final_volume=self.final_volume,
            tail_bases=self.tail_bases,
            tail_length=self.tail_length,
            repair_3prime=0.0,
        )
        pool = ligation.process(pool)

        # second step: bead-based cleanup
        cleanup = BeadCleanup(
            final_volume=self.final_volume,
        )
        pool = cleanup.process(pool)
        return pool
    
