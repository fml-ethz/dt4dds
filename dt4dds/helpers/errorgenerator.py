import numpy as np
import random
from numba import jit

from . import tools
from . import config
from .. import datastructures




@jit(nopython=True)
def arrange_and_order_with_other_list(index_list, events_list, others_list):
    result = []
    ix = 0
    for i in range(len(events_list)):
        event = events_list[i]
        indexes = index_list[ix:ix+event]
        others = others_list[ix:ix+event]
        ix += event

        result.append(sorted(zip(indexes, others), reverse=True))
    return result


@jit(nopython=True)
def arrange_and_order(index_list, events_list):
    result = []
    ix = 0
    for i in range(len(events_list)):
        event = events_list[i]
        indexes = index_list[ix:ix+event]
        ix += event

        result.append(sorted(indexes, reverse=True))
    return result


@jit(nopython=True)
def arrange_and_order_substitutions(positions_0, positions_1, positions_2, positions_3, bases_0, bases_1, bases_2, bases_3, events):
    ix = [0, 0, 0, 0]
    positions = [positions_0, positions_1, positions_2, positions_3]
    bases = [bases_0, bases_1, bases_2, bases_3]
    results = []
    for jx in range(events.shape[0]):
        event = events[jx, :]
        i_result = []
        for kx in range(4):
            i_event = event[kx]
            if not i_event: continue
            i_position = positions[kx][ix[kx]:ix[kx]+i_event]
            i_base = bases[kx][ix[kx]:ix[kx]+i_event]
            ix[kx] += i_event
            i_result.extend(list(zip(i_position, i_base)))
        results.append(sorted(i_result, reverse=True))

    return results









class ErrorGenerator():

    BLOCK_SIZE = 25000

    def __init__(self, **kwargs):
        
        length_bias = kwargs.get('length_bias', None)
        length_bias_params = kwargs.get('length_bias_params', None)

        if length_bias == 'linear':
            self.use_length_bias = True
            self.length_bias = self._linear_length_bias(length_bias_params)
        elif length_bias == 'exponential':
            self.use_length_bias = True
            self.length_bias = self._exponential_length_bias(length_bias_params)
        elif length_bias == 'points':
            self.use_length_bias = True
            self.length_bias = np.array(length_bias_params)
        else:
            self.use_length_bias = False
            self.length_bias = None


    def _linear_length_bias(self, params):
        probs = params['slope']*(np.arange(params['max_length'])+1) + params['offset']
        return probs

    def _exponential_length_bias(self, params):
        probs = np.exp(params['slope']*(np.arange(params['max_length'])+1)) + params['offset']
        return probs



    def __call__(self, sequence: datastructures.Seq, count: int):

        if not self.mean_rate > 0:
            return datastructures.SeqPool(pool_data = {str(sequence): count})

        oligo_pool = datastructures.SeqPool()

        p_ge_1error = 1-(1-self.mean_rate)**len(sequence)
        if p_ge_1error == 0:
            threshold = max(config.error_generation_max_coverage_threshold, 2*config.error_generation_error_coverage)
        else:
            threshold = config.error_generation_error_coverage / p_ge_1error
            threshold = max(
                min(threshold, config.error_generation_max_coverage_threshold), 
                2*config.error_generation_error_coverage
            )

        # limit the error generation to the specified threshold
        scale_factor = max(int(count/threshold), 1)
        if scale_factor > 1:
            offset = count % scale_factor
            oligo_pool.add_sequence(sequence, offset)
            count = (count-offset) // scale_factor

        # 
        if count > self.BLOCK_SIZE:
            count_blocks = [self.BLOCK_SIZE]*(count//self.BLOCK_SIZE)
            count_blocks.append(count % self.BLOCK_SIZE)
        else:
            count_blocks = [count]

        # 
        for count_i in count_blocks:

            events, any_events_ix = self._generate_events(sequence, count_i)

            oligo_pool.add_sequence(sequence, scale_factor*(count_i-len(any_events_ix)))

            # terminate if no events generated
            if len(any_events_ix) == 0: 
                continue

            events = events[any_events_ix]
            sum_events = np.sum(events, axis=0)

            oligos = self._generate_oligos(sequence, events, sum_events)
            # add a single copy of the mutated oligo
            oligo_pool.add_single_sequences(oligos, count=scale_factor)

        return oligo_pool




















class Substitutions(ErrorGenerator):

    OTHER_BASES_IX = [
        [1, 2, 3],
        [0, 2, 3],
        [0, 1, 3],
        [0, 1, 2]
    ]

    def __init__(self, mean_rate, bias, **kwargs):
        super().__init__(**kwargs)

        self.mean_rate = mean_rate
        self.bias = bias

        self.per_base_rate = np.array([4*mean_rate*sum(bias[base].values()) for base in tools.BASES])

        self.substituting_bias_by_base = {
            old_base: [
                bias[old_base][new_base]/sum(bias[old_base].values()) if sum(bias[old_base].values()) > 0 else 0
                for new_base in tools.base_choices(old_base)
            ] for old_base in tools.BASES
        }


    def _generate_events(self, sequence: datastructures.Seq, count: int):
        base_composition = sequence.base_composition(list=True, U2T=True)
        events_per_base = tools.rng.binomial(base_composition, self.per_base_rate, size=(count, 4))
        any_events_ix = np.where(np.any(events_per_base, axis=1))[0]

        return events_per_base, any_events_ix


    def _generate_oligos(self, sequence: datastructures.Seq, events_per_base, sum_events_per_base):

        base_idx = sequence.base_indexes(list=True, U2T=True)

        # there is no standard positional bias
        bias = [np.ones(len(idx), dtype="float") for idx in base_idx]

        # if needed, add the length bias
        if self.use_length_bias:
            bias = [bias[i]*self.length_bias[idx] for i, idx in enumerate(base_idx)]

        # determine which of the bases is substituted
        base_positions = [np.array(random.choices(
            base_idx[i],
            weights=bias[i],
            k=sum_events_per_base[i]
        ), dtype="int32") if sum_events_per_base[i] > 0 else np.array([], dtype="int32") for i in range(4)]

        # decide on the new base based on the bias
        new_base = [np.array(random.choices(
            self.OTHER_BASES_IX[i], 
            weights=self.substituting_bias_by_base[base], 
            k=sum_events_per_base[i]
        ), dtype="int32") if sum_events_per_base[i] > 0 else np.array([], dtype="int32") for i, base in enumerate(tools.BASES)]

        map_list = arrange_and_order_substitutions(*base_positions, *new_base, events_per_base)
        return map(sequence.substitution, map_list)

        








class Insertions(ErrorGenerator):

    def __init__(self, mean_rate, bias, **kwargs):
        super().__init__(**kwargs)

        self.mean_rate = mean_rate
        self.bias = np.array([bias[base] for base in tools.BASES])


    def _generate_events(self, sequence: datastructures.Seq, count: int):
        insertion_events = tools.rng.binomial(len(sequence), self.mean_rate, size=count)
        any_events_ix = np.where(insertion_events)[0]

        return insertion_events, any_events_ix



    def _generate_oligos(self, sequence: datastructures.Seq, events_per_base, sum_events_per_base):

        # there is no standard positional bias
        total_bias = np.ones(len(sequence), dtype="float")

        # if needed, add the length bias
        if self.use_length_bias:
            total_bias *= self.length_bias[np.arange(len(sequence))]

        # bias must sum to 1
        total_bias /= np.sum(total_bias)

        # for all generated events, choose a base position to insert to
        base_positions = np.array(random.choices(
            range(len(sequence)), 
            k=sum_events_per_base,
            weights=total_bias,
        ), dtype="int32")

        # 
        new_bases = np.array(random.choices(
            range(4), 
            k=sum_events_per_base, 
            weights=self.bias
        ), dtype="int32")

        map_list = arrange_and_order_with_other_list(base_positions, events_per_base, new_bases)
        return map(sequence.insertion, map_list)






class Deletions(ErrorGenerator):


    def __init__(self, mean_rate, bias, **kwargs):
        super().__init__(**kwargs)

        self.mean_rate = mean_rate
        self.mean_per_base_rate = np.array([4*mean_rate*bias[base] for base in tools.BASES])
        self.bias = bias



    def _generate_events(self, sequence: datastructures.Seq, count: int):

        # find the average per base rate for this oligo
        base_composition = sequence.base_composition(list=True, U2T=True)
        mean_oligo_rate = sum(base_composition*self.mean_per_base_rate/sum(base_composition))

        # generate random distribution of the amount of events per oligo
        events_per_base = tools.rng.binomial(len(sequence), mean_oligo_rate, size=count)

        any_events_ix = np.where(events_per_base)[0]

        return events_per_base, any_events_ix



    def _generate_oligos(self, sequence: datastructures.Seq, events_per_base, sum_events_per_base):

        # position selection is biased from the base identity
        total_bias = np.array([self.bias[b] for b in sequence])

        # if needed, add the length bias
        if self.use_length_bias:
            total_bias *= self.length_bias[np.arange(len(sequence))]

        # bias must sum to 1
        total_bias /= np.sum(total_bias)

        # for all generated events, choose a base position to delete
        base_positions = np.array(random.choices(
            range(len(sequence)), 
            k=sum_events_per_base,
            weights=total_bias,
        ), dtype="int32")

        map_list = arrange_and_order(base_positions, events_per_base)
        return map(sequence.deletion, map_list)





class DeletionEvents(ErrorGenerator):


    def __init__(self, mean_rate, bias, **kwargs):
        super().__init__(**kwargs)

        self.mean_rate = mean_rate
        self.mean_per_base_rate = np.array([4*mean_rate*bias[base] for base in tools.BASES])
        self.bias = bias


        read_bias = kwargs.get('read_bias', None)
        read_bias_params = kwargs.get('read_bias_params', None)
        if read_bias == 'points':
            self.use_read_bias = True
            self.read_bias = np.array(read_bias_params)/np.sum(read_bias_params)
        else:
            self.use_read_bias = False
            self.read_bias = None

        repeat_bias = kwargs.get('repeat_bias', None)
        repeat_bias_params = kwargs.get('repeat_bias_params', None)
        if repeat_bias == 'points':
            self.use_repeat_bias = True
            self.repeat_bias = np.array(repeat_bias_params)/np.sum(repeat_bias_params)
        else:
            self.use_repeat_bias = False
            self.repeat_bias = None
        


    def _generate_events(self, sequence: datastructures.Seq, count: int):

        # find the average per base rate for this oligo
        base_composition = sequence.base_composition(list=True, U2T=True)
        mean_oligo_rate = sum(base_composition*self.mean_per_base_rate/sum(base_composition))
        if self.use_length_bias:
            mean_oligo_rate *= self.length_bias[0:len(sequence)].mean()

        # generate random distribution of the amount of events per oligo
        if self.use_read_bias:
            events_per_oligo = tools.rng.choice(len(self.read_bias), replace=True, p=self.read_bias, size=count)
        else:
            events_per_oligo = tools.rng.binomial(len(sequence), mean_oligo_rate, size=count)

        # select only oligos with at least one event
        any_events_ix = np.where(events_per_oligo)[0]

        return events_per_oligo, any_events_ix



    def _generate_oligos(self, sequence: datastructures.Seq, events_per_base, sum_events_per_base):

        # position selection is biased from the base identity
        total_bias = np.array([self.bias[b] for b in sequence])

        # if needed, add the length bias
        if self.use_length_bias:
            total_bias *= self.length_bias[0:len(sequence)]

        # bias must sum to 1
        total_bias /= np.sum(total_bias)

        # for all generated events, choose a base position to delete
        new_positions = np.array(random.choices(
            range(len(sequence)), 
            k=sum_events_per_base,
            weights=total_bias,
        ), dtype="int32")

        # determine length of each delevent
        if self.use_repeat_bias:
            event_lengths = np.array(random.choices(
                range(len(self.repeat_bias)),
                weights=self.repeat_bias,
                k=sum_events_per_base
            ), dtype="int32")
        else:
            event_lengths = np.zeros(sum_events_per_base, dtype="int32")

        map_list = arrange_and_order_with_other_list(new_positions, events_per_base, event_lengths)
        return map(sequence.delevent, map_list)
