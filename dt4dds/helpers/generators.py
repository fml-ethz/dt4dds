from .. import datastructures
from . import tools



def create_random_sequences(n=100, l=100):
    return [''.join(tools.rng.choice(tuple(tools.BASES), size=l)) for _ in range(n)]


def attach_primers_to_sequences(seq_list, forward_primer, backward_primer):
    forward_primer_fw = datastructures.Seq(forward_primer)
    backward_primer_rc = datastructures.Seq(backward_primer).reverse_complement()

    return [''.join(map(str, [forward_primer_fw, seq, backward_primer_rc])) for seq in seq_list]


def sequences_to_fixed_mass_pool(seqlist, mass=0.5):
    pool = datastructures.SeqPool()
    pool.add_single_sequences(seqlist, count=1)
    counts = int(mass/pool.mass)-1
    pool.add_single_sequences(seqlist, count=counts)
    return pool




def create_random_pool(n=100, l=100):
    pool = datastructures.SeqPool()
    seqs = [''.join(tools.rng.choice(tuple(tools.BASES), size=l)) for _ in range(n)]
    pool.add_single_sequences(seqs, count=1)
    return pool



def create_random_pool_with_random_abundance(n=100, l=100, c_lim=[1, 100]):
    pool = datastructures.SeqPool()
    seqs = [''.join(tools.rng.choice(tuple(tools.BASES), size=l)) for _ in range(n)]
    counts = tools.rng.choice(range(c_lim[0], c_lim[1]+1), size=n)
    pool.add_sequences(seqs, counts)
    return pool


def create_random_pool_with_fixed_mass(n=100, l=100, mass=0.5):
    pool = datastructures.SeqPool()
    seqs = [''.join(tools.rng.choice(tuple(tools.BASES), size=l)) for _ in range(n)]
    pool.add_single_sequences(seqs, count=1)
    counts = int(mass/pool.mass)-1
    pool.add_single_sequences(seqs, count=counts)
    return pool




def attach_primers_to_pool(pool, forward_primer, backward_primer):
    new_pool = datastructures.SeqPool(is_doublestranded=pool.is_doublestranded, volume=pool.volume)

    forward_primer_fw = datastructures.Seq(forward_primer)
    backward_primer_rc = datastructures.Seq(backward_primer).reverse_complement()

    for sequence, count in pool:
        new_seq = "".join(map(str, [forward_primer_fw, sequence, backward_primer_rc]))
        new_pool.add_sequence(new_seq, count)

    return new_pool




def create_pool_from_sequences(seqlist, n=100):
    pool = datastructures.SeqPool()
    pool.add_single_sequences(seqlist, count=n)
    return pool