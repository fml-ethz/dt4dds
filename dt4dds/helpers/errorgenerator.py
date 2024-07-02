import numpy as np
import random
import numba
import itertools
import tqdm
import dataclasses
import time

from . import tools
from . import config
from .. import datastructures
from .. import properties

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


@numba.jit(nopython=True, cache=True)
def _extend_positions_for_events(events, positions, lengths, max_pos):
    # we extend the positions array by the number of additional consecutive errors
    new_positions = []

    ix = 0
    for i, n_events in enumerate(events):
        # we increase the number of events for each oligo by the number of additional consecutive errors
        length_of_events = lengths[ix:ix+n_events]

        pos = set()
        # we extend the positions array by the positions of these additional consecutive errors
        for j, l in enumerate(length_of_events):
            p = positions[ix+j] + np.arange(0, l+1)
            pos.update(set(p[p < max_pos]))
        
        new_positions.extend(pos)
        events[i] = len(pos)
        ix += n_events

    return events, np.array(new_positions), events.sum()


alph_as_num = np.array(list('ACGT')).view(np.int32)
base2ix = np.zeros((alph_as_num.max()+1), dtype=np.uint8)
base2ix[alph_as_num] = np.arange(4)
def seq2ix(sequence):
    """  """
    indices = base2ix[np.array([sequence,]).view(np.int32)]
    return indices


ix2base = np.array(['A', 'C', 'G', 'T'], dtype='<U1')
def ix2seq(list_of_indices, equal_length=False):
    """  """
    if equal_length:
        return ix2base[list_of_indices].view(f'U{list_of_indices[0].size}').flatten()
    else:
        return [ix2base[i].view(f'U{i.size}')[0] for i in list_of_indices]






@dataclasses.dataclass
class ErrorGenerator():

    # user-determined error parameters
    mean_rate: float = 0
    bias: dict = None
    length_bias: list = None
    read_bias: list = None
    repeat_bias: list = None

    # parameters for coverage of generated erroneous sequences
    error_generation_error_coverage: int = config.error_generation_error_coverage
    error_generation_max_coverage_threshold: int = config.error_generation_max_coverage_threshold

    # whether mutated oligos have length identical to parent, speeds up sequence generation
    mutated_oligos_have_equal_length = False

    # name of the error generator
    name = "ErrorGenerator"

    def __repr__(self):
        return f"ErrorGenerator({self.name})"

    # 
    # initializations
    # 

    def _initialize_oligo_biases(self):
        ''' Initializes the oligo biases. '''

        # length bias gives error distribution between positions in oligos
        if self.length_bias is not None:
            self.length_bias = np.array(self.length_bias)/np.sum(self.length_bias)
            self.use_length_bias = True
            logger.debug(f"{self}: Using length bias for error generation.")
        else:
            self.use_length_bias = False

        # read bias gives error distribution between individual reads
        if self.read_bias is not None:
            self.read_bias = np.array(self.read_bias)/np.sum(self.read_bias)
            self.use_read_bias = True
            logger.debug(f"{self}: Using read bias for error generation. This overrides the mean error rate for error generation.")
        else:
            self.use_read_bias = False

        # repeat bias gives error repeat distribution for consecutive events
        if self.repeat_bias is not None:
            self.repeat_bias = np.array(self.repeat_bias)/np.sum(self.repeat_bias)
            self.use_repeat_bias = True
            logger.debug(f"{self}: Using repeat bias for error generation.")
        else:
            self.use_repeat_bias = False


    def _initialize_repeat_iterator(self):
        ''' Initializes the iterators for the event lengths. '''
        event_lengths = random.choices(
            np.arange(len(self.repeat_bias)),
            k=1000,
            weights=self.repeat_bias,
        )
        self.event_lengths = itertools.cycle(event_lengths)


    def _initialize_read_iterator(self):
        ''' Initializes the iterators for the read event count. '''
        event_reads = random.choices(
            np.arange(len(self.read_bias)),
            k=1000,
            weights=self.read_bias,
        )
        self.event_reads = itertools.cycle(event_reads)


    # 
    # iterator getters
    # 

    def _get_event_lengths(self, n_events: int):
        ''' Generates the number of consecutive events for n_events events. '''
        return np.fromiter(self.event_lengths, dtype="H", count=n_events)
    

    def _get_event_counts(self, n_oligos: int):
        ''' Generates the count of events for n_oligos oligos. '''
        return np.fromiter(self.event_reads, dtype="H", count=n_oligos)
    

    # 
    # event generation
    # 

    def _generate_events(self, sequence: np.array, count: int):
        ''' Generates the number of error events per oligo for mutation. '''

        # generate random distribution of the amount of events per oligo
        if self.use_read_bias:
            events_per_oligo = self._get_event_counts(count)
        else:
            # correct for the length bias if required
            if self.use_length_bias:
                mean_oligo_rate = self.mean_rate * self.length_bias[0:len(sequence)].mean()
            else:
                mean_oligo_rate = self.mean_rate
            events_per_oligo = tools.rng.binomial(len(sequence), mean_oligo_rate, size=count)

        # return the generated events
        return events_per_oligo


    def _get_sequence_position_bias(self, sequence: np.array):
        
        # get the intial position bias caused by base positions
        position_bias = self.base_bias[sequence]
        
        # if needed, add the length bias
        if self.use_length_bias:
            position_bias *= self.length_bias[0:len(sequence)]

        return position_bias


    # 
    # public entry
    # 
                                    
    def mutate_pool(self, pool: datastructures.SeqPool, disable_progressbar=False):
        ''' Mutates a SeqPool. '''

        # initialize the biases
        self._initialize_oligo_biases()
        self._initialize_base_biases()

        # short-circuit if no errors expected
        if not self.mean_rate > 0:
            logger.debug(f"{self}: Skipping, mean error rate is zero.")
            return pool
        if self.use_read_bias and self.read_bias[0] == 1:
            logger.debug(f"{self}: Skipping, all oligos are error-free.")
            return pool
        
        # time the mutation process
        t = time.time()
        
        # initialize the iterators
        self._initialize_base_iterator()
        if self.use_repeat_bias: self._initialize_repeat_iterator()
        if self.use_read_bias: self._initialize_read_iterator()
                
        # generate the mutated pool
        mutated_pool = datastructures.SeqPool()
        pool_generator = tqdm.tqdm(pool, desc=f"Generating {self.name}", disable=disable_progressbar or not config.show_progressbars, total=pool.n_sequences)
        for seq, count in pool_generator:
            self._add_mutations(seq.sequence, count, mutated_pool)

        logger.info(f"{self} finished, taking {time.time()-t:.2f} seconds for {pool.n_oligos} oligos.")
        return mutated_pool


    def mutate_iterator(self, iterator: list, disable_progressbar=False):
        ''' Mutates a iterator yielding SeqPool objects. '''

        # initialize the biases
        self._initialize_oligo_biases()
        self._initialize_base_biases()

        # short-circuit if no errors expected
        if not self.mean_rate > 0:
            logger.debug(f"{self}: Skipping, mean error rate is zero.")
            for pool in iterator: yield pool
            return
        if self.use_read_bias and self.read_bias[0] == 1:
            logger.debug(f"{self}: Skipping, all oligos are error-free.")
            for pool in iterator: yield pool
            return
        
        # time the mutation process
        t = time.time()
        
        # initialize the iterators
        self._initialize_base_iterator()
        if self.use_repeat_bias: self._initialize_repeat_iterator()
        if self.use_read_bias: self._initialize_read_iterator()

        # generate the mutated pools
        pool_generator = tqdm.tqdm(iterator, desc=f"Generating {self.name}", disable=disable_progressbar or not config.show_progressbars)
        tot_n_oligos = 0
        for pool in pool_generator:
            mutated_pool = datastructures.SeqPool()
            for seq, count in pool:
                self._add_mutations(seq.sequence, count, mutated_pool)
                tot_n_oligos += count
            yield mutated_pool
        logger.info(f"{self} finished, taking {time.time()-t:.2f} seconds for {tot_n_oligos} oligos.")


    # 
    # related to mutation generation
    # 

    def _get_positions(self, sequence: datastructures.Seq, events_per_oligo: np.array, sum_events_per_oligo: int):

        # position selection is biased from the base identity
        position_bias = self._get_sequence_position_bias(sequence)

        # for all generated events, choose base positions for the errors
        base_positions = np.array(random.choices(
            range(len(sequence)), 
            k=sum_events_per_oligo,
            weights=position_bias,
        ), dtype="H")

        # if we have repeat events, we need to adjust the positions and event counts
        if self.use_repeat_bias:
            # choose the length of the insertion based on the repeat bias
            event_lengths = self._get_event_lengths(sum_events_per_oligo)

            # correct the positions and events for the additional consecutive errors
            events_per_oligo, base_positions, sum_events_per_oligo = _extend_positions_for_events(
                events_per_oligo, 
                base_positions, 
                event_lengths,
                len(sequence)
            )

        return events_per_oligo, base_positions, sum_events_per_oligo
        

    def _get_coverage(self, seq_len: int):
        ''' Generates the error coverage for the number of oligos to generate. '''
        
        # get the expected fraction of oligos with at least one error
        if self.use_read_bias:
            # probability of at least one error is (1 - prob. of no error)
            p_ge_1error = 1 - self.read_bias[0]
        else:
            # probability of at least one error is (1 - prob. of no error)
            p_ge_1error = 1-(1-self.mean_rate)**seq_len

        # to get sufficient oligos with at least one error, we need to generate at least this many oligos
        coverage = int(self.error_generation_error_coverage / p_ge_1error) + 1

        # but we don't want to generate more than the specified maximum threshold
        return min(coverage, self.error_generation_max_coverage_threshold), p_ge_1error
    

    def _add_mutations(self, sequence: str, count: int, mutated_pool: datastructures.SeqPool):
        ''' Generates the number of error events per oligo for mutation. '''

        # get the coverage for the number of oligos to generate
        coverage, p_ge_1error = self._get_coverage(len(sequence))

        # establish an integer scaling factor to limit the error generation to the coverage
        scale_factor = max(count // coverage, 1)
        scale_count = min(count, coverage)
        scale_overflow = count - scale_count*scale_factor
        scale_overflow_w_errors = int(scale_overflow*p_ge_1error)
        scale_overflow_wo_errors = scale_overflow - scale_overflow_w_errors

        # convert sequence into an array of indices
        ix_sequence = seq2ix(sequence)

        # generate error events
        events = self._generate_events(ix_sequence, scale_count)

        # get the indices and amount of oligos with at least one error
        any_events_ix = np.where(events)[0]
        n_oligos_with_events = len(any_events_ix)

        # add the unmutated oligos
        mutated_pool.add_sequence(sequence, scale_factor*(scale_count-n_oligos_with_events) + scale_overflow_wo_errors)

        # terminate if no events occured
        if n_oligos_with_events == 0: 
            mutated_pool.add_sequence(sequence, scale_overflow_w_errors)
            return
        
        # get the entries where events occured and their total number
        events = events[any_events_ix]
        n_events = sum(events)

        # generate and add the mutated oligos
        events, positions, n_events = self._get_positions(ix_sequence, events, n_events)
        new_bases = self._get_new_bases(n_events)
        new_indices = self._mutate(ix_sequence, events, positions, new_bases, n_events)

        # convert the indices back to sequences
        oligos = ix2seq(new_indices, equal_length=self.mutated_oligos_have_equal_length)
        mutated_pool.add_single_sequences(oligos, count=scale_factor)

        # copy the efficiency property for all new oligos
        properties.AmplificationEfficiency.duplicate_efficiencies(sequence, oligos)

        # for the overflow, select a subset of the mutated oligos
        if scale_overflow_w_errors > 0:
            mutated_pool.add_single_sequences(random.choices(oligos, k=scale_overflow_w_errors), count=1)
        return

    # 
    # implemented by child classes
    # 
    
    def _initialize_base_iterator(self):
        raise NotImplementedError
    
    def _get_new_bases(self, n_bases: int):
        raise NotImplementedError

    @staticmethod
    def _mutate(self, sequence: datastructures.Seq, events_per_oligo, base_positions, sum_events_per_oligo):
        raise NotImplementedError

    def _initialize_base_biases(self):
        raise NotImplementedError



class Substitutions(ErrorGenerator):
    ''' Implements substitutions / subevents. '''

    mutated_oligos_have_equal_length = True
    name = "Substitutions"

    OTHER_BASES_IX = [
        [1, 2, 3],
        [0, 2, 3],
        [0, 1, 3],
        [0, 1, 2]
    ]

    def _initialize_base_biases(self):

        # the bias for any substitution occuring for an original base
        self.base_bias = np.array([sum(self.bias[base].values()) for base in tools.BASES])

        # for each base, we compute the bias for the substituting base
        self.substituting_base_bias = np.array([
            [
                self.bias[old_base][new_base]/sum(self.bias[old_base].values()) if sum(self.bias[old_base].values()) > 0 else 0
                for new_base in tools.base_choices(old_base)
            ] for old_base in tools.BASES
        ])
    

    def _initialize_base_iterator(self):
        ''' Initializes the iterators for the new bases. '''
        new_bases = [random.choices(
            self.OTHER_BASES_IX[i], 
            weights=self.substituting_base_bias[i], 
            k=100
        ) for i in range(4)]
        self.new_bases = [itertools.cycle(new_base) for new_base in new_bases]


    def _get_new_bases(self, n_bases: int):
        ''' Generates the new bases for n_bases bases. '''
        return np.array([np.fromiter(self.new_bases[i], dtype="H", count=n_bases) for i in range(4)], dtype="H")


    @staticmethod
    @numba.jit(nopython=True, cache=True)
    def _mutate(sequence: np.array, events_per_oligo: np.array, base_positions: np.array, new_bases: np.array, sum_events_per_oligo: int):
        '''  '''
        n_seqs = events_per_oligo.shape[0]
        new_sequences = np.ascontiguousarray(sequence.repeat(n_seqs).reshape((-1, n_seqs)).T)

        ix = 0
        for i in range(len(events_per_oligo)):
            n_events = events_per_oligo[i]
            # pull an adequate number of positions
            positions = base_positions[ix:ix+n_events]

            # for each position, grab an adequate new base and substitute
            for j in range(len(positions)):
                pos = positions[j]
                base = new_bases[sequence[pos], ix+j]
                new_sequences[i, pos] = base

            # increment the index
            ix += n_events

        # return all new oligos
        return new_sequences



class Insertions(ErrorGenerator):
    ''' Implements insertions / insevents. '''

    mutated_oligos_have_equal_length = False
    name = "Insertions"

    def _initialize_base_biases(self):
        # we assume there's no inherent base bias for insertions
        self.base_bias = np.ones(4, dtype="float")
        # for each base, we compute the bias for the inserted base
        self.inserted_base_bias = np.array([self.bias[base] for base in tools.BASES])


    def _initialize_base_iterator(self):
        ''' Initializes the iterators for the new bases. '''
        new_bases = random.choices(
            np.arange(4), 
            weights=self.inserted_base_bias, 
            k=100
        )
        self.new_bases = itertools.cycle(new_bases)


    def _get_new_bases(self, n_bases: int):
        ''' Generates the new bases for n_bases bases. '''
        return np.fromiter(self.new_bases, dtype="H", count=n_bases)


    @staticmethod
    @numba.jit(nopython=True, cache=True)
    def _mutate(sequence: np.array, events_per_oligo: np.array, base_positions: np.array, new_bases: np.array, sum_events_per_oligo: int):
        '''  '''
        new_sequences = []
        orig_len = len(sequence)

        ix = 0
        for i in range(len(events_per_oligo)):
            n_events = events_per_oligo[i]
            # pull an adequate number of positions and bases
            positions = base_positions[ix:ix+n_events]
            bases = new_bases[ix:ix+n_events]
            
            # mask the sequence array
            new_sequence = np.zeros(orig_len+n_events, dtype="H")

            # copy the original sequence
            new_sequence[:orig_len] = sequence

            # insert the new bases
            for pos, new_base in sorted(zip(positions, bases), reverse=True):
                new_sequence[pos+2:] = new_sequence[pos+1:-1]
                new_sequence[pos+1] = new_base
            
            new_sequences.append(new_sequence)

            # increment the index
            ix += n_events

        # return all new oligos
        return new_sequences



class Deletions(ErrorGenerator):
    ''' Implements deletions / delevents. '''

    mutated_oligos_have_equal_length = False
    name = "Deletions"

    def _initialize_base_biases(self):
        # we assume each base has a constant, inherent bias for deletions
        self.base_bias = np.array([self.bias[base] for base in tools.BASES])
    

    def _initialize_base_iterator(self):
        ''' Initializes the iterators for the new bases. '''
        # deletions don't have new bases
        pass


    def _get_new_bases(self, n_bases: int):
        ''' Generates the new bases for n_bases bases. '''
        # deletions don't have new bases
        return np.array([], dtype="H")
    

    @staticmethod
    @numba.jit(nopython=True, cache=True)
    def _mutate(sequence: np.array, events_per_oligo: np.array, base_positions: np.array, new_bases: np.array, sum_events_per_oligo: int):
        '''  '''
        new_sequences = []
        mask = np.ones(sequence.size, dtype=np.bool_)

        ix = 0
        for i in range(len(events_per_oligo)):
            n_events = events_per_oligo[i]
            # pull an adequate number of positions
            positions = base_positions[ix:ix+n_events]
            
            # mask the sequence array
            mask[:] = True
            mask[positions] = False
            new_sequence = sequence[mask]
            new_sequences.append(new_sequence)

            # increment the index
            ix += n_events

        # return all new oligos
        return new_sequences