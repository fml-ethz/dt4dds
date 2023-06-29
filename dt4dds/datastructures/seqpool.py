import csv
import gzip
from random import shuffle
import numpy as np
import pandas as pd
import itertools
import collections

from ..helpers import tools
from ..helpers import config

from .seq import Seq

import logging
logger = logging.getLogger(__name__)



class SeqPool():

    #
    # Initialization and generator methods
    #

    def __init__(self, pool_data=None, volume=None, is_doublestranded=True):

        # dict with "sequence": abundance
        self._seqdict = pool_data if pool_data else {}

        self.volume = volume
        """ Volume of this pool in uL. """

        self.is_doublestranded = is_doublestranded
        """ Whether this pool contains double-stranded or single-stranded DNA. """


    def __iter__(self):
        for sequence, count in self._seqdict.items():
            yield Seq(sequence), count


    def __contains__(self, sequence):
        return sequence in self._seqdict


    def __getitem__(self, sequence):
        return self._seqdict[sequence]


    def __setitem__(self, sequence, value):
        self._seqdict[sequence] = value


    def __repr__(self):
        df = self.to_dataframe()
        return "\n".join([
            f"SeqPool(data=",
            str(df.sort_values(by=["value"], axis=0, ascending=False)), 
            f"volume={self.volume}, is_ds={self.is_doublestranded}, "
            f"n_seqs={len(df.index)}, n_oligos={df['value'].sum()}, "
            f"mean={df['value'].mean():0.1f}, std={df['value'].std():0.1f}, "
            f"median={df['value'].median():0.1f}, mass={self.mass:0.4f})",
        ])


    #
    # Class functions
    #

    @classmethod
    def combine(cls, *other_pools):
        pool = cls()
        pool.combine_with(*other_pools)
        return pool

    @classmethod
    def from_csv(cls, filepath, sep=',', header=None):
        df_dict = pd.read_csv(filepath, sep=sep, header=header).to_dict(orient='list')
        instance = cls(pool_data = dict(zip(df_dict[0], df_dict[1])))
        return instance



    #
    # Convenience functions
    #


    def sequences(self):
        """Generator for the sequences contained in this SeqPool, each converted to a Seq object."""
        for sequence in self._seqdict.keys():
            yield Seq(sequence)


    def counts(self):
        """Generator for the counts contained in this SeqPool."""
        for count in self._seqdict.values():
            yield count


    def copy(self):
        """Returns a deep copy of this SeqPool."""
        new_instance = SeqPool(
            pool_data={sequence: value for sequence, value in self._seqdict.items()},
            volume=self.volume,
            is_doublestranded=self.is_doublestranded
        )
        return new_instance


    def combine_with(self, *other_pools):
        """Adds all sequences from the other SeqPool(s) to this SeqPool."""
        for other_pool in other_pools:
            self.add_sequences(other_pool.sequences(), other_pool.counts())
            if self.volume is not None and other_pool.volume is not None:
                self.volume += other_pool.volume
        return self


    def combine_with2(self, other_pools):
        """Adds all sequences from the other SeqPool(s) to this SeqPool."""
        for other_pool in other_pools:
            self.add_sequences(other_pool.sequences(), other_pool.counts())
            if self.volume is not None and other_pool.volume is not None:
                self.volume += other_pool.volume
        return self


    @property
    def n_sequences(self):
        """ Number of unique sequences in this pool. """
        return len(self._seqdict.keys())

    @property
    def n_oligos(self):
        """ Total number of oligos in this pool. """
        return sum(self._seqdict.values())


    @property
    def n_nucleotides(self):
        """ Total number of nucleotides in this pool. """
        lengths = np.array([len(s) for s in self.sequences()], dtype=float)
        counts = np.array([c for c in self.counts()], dtype=float)

        return int(np.sum(lengths*counts, dtype=float))

    
    @property
    def mass_concentration(self):
        """ Returns the DNA concentration in ng/uL. """
        if self.volume is None:
            raise AttributeError("SeqPool does not have a specified volume, cannot calculate concentration.")
        return self.mass/self.volume

    @property
    def molar_concentration(self):
        """ Returns the DNA concentration in nM. """
        if self.volume is None:
            raise AttributeError("SeqPool does not have a specified volume, cannot calculate concentration.")
        return 1e9*self.n_oligos/(tools.NA*self.volume)


    @property
    def mass(self):
        """ Returns the total DNA weight in nanogram. Formulas as used by NEBioCalculator."""
        if self.is_doublestranded:
            return (1e9/tools.NA) * (self.n_nucleotides * 617.96 + self.n_oligos * 36.04)
        else:
            return (1e9/tools.NA) * (self.n_nucleotides * 308.97 + self.n_oligos * 18.02)

    @property
    def moles(self):
        """"""
        return 1e9*self.n_oligos/tools.NA


    #
    # Cropping and Filtering
    #

    def simplify(self, fixed_coverage=None):

        if not fixed_coverage:
            fixed_coverage = config.seqpool_maximum_pool_diversity

        oligo_factor = self.n_oligos//fixed_coverage

        if oligo_factor > 1:
            logger.info(f"Pool simplified by a factor of {oligo_factor}")
            self._seqdict = self.sample_by_factor(1/oligo_factor, accurate=True)._seqdict
            self._seqdict.update((seq, int(count*oligo_factor)) for seq, count in self._seqdict.items())


    def split_by_threshold(self, threshold):
        """"""
        seqdict_in = {sequence: value for sequence, value in self._seqdict.items() if value >= threshold}
        pool_in = SeqPool(pool_data=seqdict_in, is_doublestranded=self.is_doublestranded)

        pool_out = self.copy()
        pool_out.remove_sequences(pool_in.sequences())
        
        return pool_in, pool_out


    def filter_by_length(self, l_min, l_max, inplace=False):
        """"""
        seqdict = {seq: value for seq, value in self._seqdict.items() if (len(seq) >= l_min) and (len(seq) <= l_max)}
        if inplace: 
            self._seqdict = seqdict
            return None
        else:
            return SeqPool(pool_data=seqdict, is_doublestranded=self.is_doublestranded)



    def sample_by_factor(self, dil_factor, remove_sampled_oligos=False, accurate=False):

        if dil_factor > 1:
            raise ArithmeticError(f"Cannot sample with {dil_factor}, must be <1.")

        sample_pool = SeqPool(is_doublestranded=self.is_doublestranded)

        if dil_factor == 0:
            return sample_pool

        if dil_factor == 1:
            sample_pool = self.copy()
            if remove_sampled_oligos: self._seqdict = {}
            return sample_pool

        counts = np.array(list(self._seqdict.values()), dtype="int64")
        n_oligos = int(dil_factor*sum(counts))
        new_counts = tools.rng.binomial(counts, dil_factor)

        # ensure count is exact
        while accurate and (diff := n_oligos - np.sum(new_counts)):

            counter = collections.Counter(tools.rng.choice(
                self.n_sequences,
                p=new_counts/np.sum(new_counts),
                size=np.abs(diff)
            ))

            diff_counts = np.zeros(len(counts), dtype="int")
            diff_counts[list(counter.keys())] = list(counter.values())

            if diff < 0:
                new_counts -= diff_counts
            else:
                new_counts += diff_counts
            new_counts = np.maximum(new_counts, 0)


        sample_pool.add_sequences(self.sequences(), new_counts)

        if remove_sampled_oligos:
            self.remove_sequences(list(self.sequences()), new_counts)

        if self.volume: sample_pool.volume = dil_factor*self.volume
        return sample_pool


    def sample_by_counts(self, n_oligos, remove_sampled_oligos=False, accurate=True):

        n_oligos = int(n_oligos)
        dil_factor = n_oligos/self.n_oligos

        sample_pool = self.sample_by_factor(dil_factor, remove_sampled_oligos=remove_sampled_oligos, accurate=accurate)
        return sample_pool


    def sample_by_volume(self, volume, remove_sampled_oligos=False, accurate=False):
        """"""
        if not self.volume:
            raise AttributeError(f"Cannot sample by volume, pool has no volume set.")

        dilution_factor = volume/self.volume

        if dilution_factor > 1:
            raise ArithmeticError(f"Cannot sample with {volume} uL, only {self.volume} uL available.")

        sample_pool = self.sample_by_factor(dilution_factor, remove_sampled_oligos=remove_sampled_oligos, accurate=accurate)
        self.volume -= volume
        sample_pool.volume = volume
        return sample_pool


    def sample_by_mass(self, mass, remove_sampled_oligos=False, accurate=False):
        """"""

        dilution_factor = mass/self.mass

        if dilution_factor > 1:
            raise ArithmeticError(f"Cannot sample with {mass} ng, only {self.mass} ng available.")

        sample_pool = self.sample_by_factor(dilution_factor, remove_sampled_oligos=remove_sampled_oligos, accurate=accurate)
        if self.volume: sample_pool.volume = dilution_factor*self.volume
        return sample_pool


    def sample_by_moles(self, moles, remove_sampled_oligos=False, accurate=False):

        dilution_factor = moles/self.moles

        if dilution_factor > 1:
            raise ArithmeticError(f"Cannot sample with {moles} nmol, only {self.moles} nmol available.")

        sample_pool = self.sample_by_factor(dilution_factor, remove_sampled_oligos=remove_sampled_oligos, accurate=accurate)
        if self.volume: sample_pool.volume = dilution_factor*self.volume
        return sample_pool


    #
    # Splitting
    #
    
    def split(self, n, sequentially=False):

        n = int(n)
        
        if sequentially:
            it = iter(self._seqdict.items())
            size = round(self.n_sequences/n)
    
            pools = [SeqPool(pool_data=dict(itertools.islice(it, size))) for _ in range(n-1)]
            pools.append(SeqPool(pool_data=dict(it)))

        else:
            pools = []
            for i in range(n):
                it = iter(self._seqdict.items())
                pools.append(SeqPool(pool_data=dict(itertools.islice(it, i, self.n_sequences, n)), is_doublestranded=self.is_doublestranded))

        return pools


    #
    # Sequence manipulation
    #

    def remove_sequence(self, sequence, value=None):
        """ Removes the specified sequence with the specified value from this SeqPool. Completely removes the specified sequence from this SeqPool if value=None. """

        if value is None:
            self._seqdict.pop(sequence, None)
            return

        current_value = self._seqdict.get(sequence, 0)
        if value > self._seqdict.get(sequence, 0):
            logging.warning(f"Attempting to remove {value} copies, but only {current_value} are present. Sequence '{sequence}'.")
        
        if current_value-value <= 0:
            self._seqdict.pop(sequence, None)
        else:
            self._seqdict[sequence] = int(current_value-value)


    def remove_sequences(self, sequences, values=None):
        if values is None: 
            for sequence in sequences:
                self.remove_sequence(sequence)
        else:
            for sequence, value in zip(sequences, values):
                self.remove_sequence(sequence, value)


    def add_sequence(self, sequence, value):
        """Adds the specified sequence with the specified value to the SeqPool."""
        if value == 0: return
        seq = str(sequence)
        self._seqdict[seq] = self._seqdict.get(seq, 0) + int(value)


    def add_sequences(self, sequences, values):
        """Adds the specified sequences with the specified values to the SeqPool."""
        sequences = map(str, sequences)
        for sequence, value in zip(sequences, values):
            if value == 0 or sequence == "": continue
            self._seqdict[sequence] = self._seqdict.get(sequence, 0) + int(value)


    def add_single_sequences(self, sequences, count=1):
        """Adds the specified sequences with the same, specified count to the SeqPool."""
        sequences = map(str, sequences)
        count = int(count)
        for seq in sequences:
            self._seqdict[seq] = self._seqdict.get(seq, 0) + count
        


    #
    # Value manipulation
    #

    def rebase_values(self, total=1.0):
        """Rebases each sequence's value so that the sum of all sequences' values are equal to the given total."""
        factor = total/sum(self._seqdict.values())
        for sequence in self._seqdict.keys():
            self._seqdict[sequence] *= factor


    #
    # Saving and exporting
    #

    def save_as_csv(self, filepath):
        """Saves this SeqPool's sequence and values in a csv file at filepath."""
        with open(filepath, 'w', newline='', encoding='utf-8') as csv_file:  
            writer = csv.writer(csv_file)
            writer.writerows([sequence, value] for sequence, value in self._seqdict.items())


    def save_as_fasta(self, filepath, subsample=None):
        """"""

        if subsample:
            seqpool = self.sample_by_counts(int(subsample))
        else:
            seqpool = self

        with open(filepath, 'w', encoding='utf-8') as f:
            for index, (sequence, count) in enumerate(seqpool):
                fasta_string = '\n'.join([
                    f'> #{str(index).zfill(6)}, abundance: {str(count)}',
                    str(sequence),
                    '\n'
                ])
                f.writelines(fasta_string)


    def save_as_fastq(self, filepath, subsample=None, paired_read=None):
        """"""

        if subsample:
            seqpool = self.sample_by_counts(int(subsample))
        else:
            seqpool = self

        f = gzip.open(filepath, 'wt')
        for index, (sequence, count) in enumerate(seqpool): 
            for i in range(count):
                string = '\n'.join([
                    f"@Seq{str(index).zfill(9)}:{str(i+1).zfill(9)}",
                    str(sequence),
                    "+",
                    "F"*len(sequence),
                    '' # final newline
                ])

                f.writelines(string)
        f.close()

        if paired_read:
            f = gzip.open(paired_read, 'wt')
            for index, (sequence, count) in enumerate(seqpool): 
                for i in range(count):
                    string = '\n'.join([
                        f"@Seq{str(index).zfill(9)}:{str(i+1).zfill(9)}",
                        str(sequence.reverse_complement()),
                        "+",
                        "F"*len(sequence),
                        '' # final newline
                    ])

                    f.writelines(string)
            f.close()


    def to_dataframe(self):
        """"""
        return pd.DataFrame.from_dict({"sequence": self._seqdict.keys(), "value": self._seqdict.values()})


    def to_list(self):
        """"""
        return [sequence for sequence, count in self for _ in range(count)]


    def to_dict(self):
        """"""
        return self._seqdict.copy()