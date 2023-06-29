import numpy as np

from ..helpers import tools

import logging
logger = logging.getLogger(__name__)


class Seq():

    __slots__ = ("sequence")


    def __init__(self, seq, normalized=False) -> None:
        self.sequence = seq if normalized else tools.standardize_seq(seq)


    def __len__(self):
        return len(self.sequence)

    
    def __hash__(self):
        return hash(self.sequence)


    def __eq__(self, other):
        return self.sequence == str(other)


    def __contains__(self, item):
        return item in self.sequence


    def __repr__(self):
        return f"Seq({self.sequence})"


    def __str__(self):
        return self.sequence


    def __iter__(self):
        return self.sequence.__iter__()


    def __getitem__(self, key):
        """Return a subsequence as a single letter or as a sequence object."""
        if isinstance(key, int):
            # Return single base
            return self.sequence[key]
        else:
            # Return subsequence as Seq instance
            return Seq(self.sequence[key], normalized=True)


    def __add__(self, other):
        return Seq(self.sequence + str(other), normalized=True)


    def __radd__(self, other):
        return Seq(str(other) + self.sequence, normalized=True)


    def copy(self):
        """ Returns a copy of this sequence. """
        return Seq(self.sequence, normalized=True)


    def replace(self, from_base, to_base):
        self.sequence = self.sequence.replace(from_base, to_base)
    

    def substitution(self, events):
        """"""
        sequence = self.sequence
        for index, new_base in events:
            sequence = sequence[0:index] + tools.BASES[new_base] + sequence[index+1:]
        return Seq(sequence, normalized=True)


    def insertion(self, events):
        """"""
        sequence = self.sequence
        for index, new_base in events:
            sequence = sequence[0:index+1] + tools.BASES[new_base] + sequence[index+1:]
        return Seq(sequence, normalized=True)


    def deletion(self, indexes):
        """"""
        sequence = self.sequence
        for index in indexes:
            sequence = sequence[0:index] + sequence[index+1:]
        return Seq(sequence, normalized=True)


    def delevent(self, events):
        """"""
        sequence = self.sequence
        for index, length in events:
            sequence = sequence[0:index] + sequence[index+length+1:]
        return Seq(sequence, normalized=True)


    def find(self, subseq):
        """ Gives the index in this sequence at which the subseq starts, or -1 if no match. """
        return self.sequence.find(str(subseq))


    def find_primer(self, primerseq, min_overlap=10):
        """ Identifies the index in this sequence at which the last min_overlap bases of the primer match. Returns None if no match found. """
        if len(primerseq) < min_overlap:
            logger.warning(f"The used primer {primerseq} is shorter than the minimum primer length.")
            return None
        subseq = str(primerseq)[-min_overlap:]
        ix = self.sequence.find(subseq)
        if ix >= 0:
            return ix + min_overlap
        else:
            return None


    def find_primers(self, primers_list, min_overlap=10):
        """  """
        for primerseq in primers_list:
            if (index_start := self.find_primer(primerseq, min_overlap=min_overlap)):
                break
        if index_start:
            return index_start, primerseq
        else:
            return None, None


    def reverse_complement(self):
        """Return the reverse complement as a new sequence."""
        seq_rc = ''.join(tools.BASE_COMPLEMENT[c] for c in reversed(self.sequence))
        return Seq(seq_rc, normalized=True)


    def reverse(self):
        """Return the reverse as a new sequence."""
        seq_r = self.sequence[::-1]
        return Seq(seq_r, normalized=True)


    def complement(self):
        """Return the complement as a new sequence."""
        seq_c = ''.join(tools.BASE_COMPLEMENT[c] for c in self.sequence)
        return Seq(seq_c, normalized=True)


    def base_composition(self, list=False, U2T=False):
        """  """
        if U2T:
            if list:
                result = np.array([self.sequence.count(base) for base in tools.BASES])
                result[tools.BASE2IDX["T"]] += self.sequence.count("U")
                return result
            else:
                result = {base: self.sequence.count(base) for base in tools.BASES}
                result["T"] += self.sequence.count("U")
                return result

        else:
            if list:
                return np.array([self.sequence.count(base) for base in tools.BASES_WITH_U])
            else:
                return {base: self.sequence.count(base) for base in tools.BASES_WITH_U}


    def base_indexes(self, list=False, U2T=False):
        """  """
        if U2T:
            indexes = [[] for _ in range(4)]
            for i, base in enumerate(self.sequence):
                indexes[tools.BASE2IDX_U2T[base]].append(i)
            if not list: return {base: indexes[i] for i, base in enumerate(tools.BASES)}
            return indexes

        else:
            indexes = [[] for _ in range(5)]
            for i, base in enumerate(self.sequence):
                indexes[tools.BASE2IDX[base]].append(i)
            if not list: return {base: indexes[i] for i, base in enumerate(tools.BASES_WITH_U)}
            return indexes