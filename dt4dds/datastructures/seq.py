from ..helpers import tools

import logging
logger = logging.getLogger(__name__)

# https://bioinformatics.stackexchange.com/questions/3583/what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho
tab = str.maketrans("ACTG", "TGAC")


class Seq():

    __slots__ = ("sequence", "damaged_3prime")

    def __init__(self, seq) -> None:
        if seq[-1] == "_":
            seq = seq[:-1]
            self.damaged_3prime = True
        else:
            self.damaged_3prime = False
        self.sequence = str(seq)


    def __len__(self):
        return len(self.sequence)

    
    def __hash__(self):
        return hash(self.sequence) + int(self.damaged_3prime)


    def __eq__(self, other):
        return str(self) == str(other)


    def __contains__(self, item):
        return item in self.sequence


    def __repr__(self):
        return f"Seq({self.sequence})"


    def __str__(self):
        if self.damaged_3prime:
            return str(self.sequence) + "_"
        else:
            return str(self.sequence)


    def __iter__(self):
        return self.sequence.__iter__()


    def __getitem__(self, key):
        """Return a subsequence as a single letter or as a sequence object."""
        if isinstance(key, int):
            # Return single base
            return self.sequence[key]
        else:
            # Return subsequence as Seq instance
            return Seq(self.sequence[key])


    def __add__(self, other):
        return Seq(self.sequence + str(other))


    def __radd__(self, other):
        return Seq(str(other) + self.sequence)


    def copy(self):
        """ Returns a copy of this sequence. """
        return Seq(str(self))
    

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
        
    def rfind_primer(self, primerseq, min_overlap=10):
        """ Identifies the index in this sequence at which the last min_overlap bases of the primer match. Returns None if no match found. """
        if len(primerseq) < min_overlap:
            logger.warning(f"The used primer {primerseq} is shorter than the minimum primer length.")
            return None
        subseq = str(primerseq)[:min_overlap]
        ix = self.sequence.rfind(subseq)
        if ix >= 0:
            return ix
        else:
            return None


    def reverse_complement(self):
        """Return the reverse complement as a new sequence."""
        return Seq(self.sequence.translate(tab)[::-1])


    def reverse(self):
        """Return the reverse as a new sequence."""
        return Seq(self.sequence[::-1])


    def complement(self):
        """Return the complement as a new sequence."""
        return Seq(self.sequence.translate(tab))