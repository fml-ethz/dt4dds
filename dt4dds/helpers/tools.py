import re
import numpy as np
import pathlib
import gzip
import functools
import Bio.SeqIO


# random instance for all random operations, change via tools.set_seed(seed)
rng = np.random.default_rng(np.random.SFC64())

BASES = ["A", "C", "G", "T"]
BASES_WITH_U = ["A", "C", "G", "T", "U"]
BASE_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}
BASE2IDX = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 4}
BASE2IDX_U2T = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}

BASES_MW = {"A": 313.2, "C": 289.2, "G": 329.2, "T": 304.2, "U": 290.2} # g/mol

NA = 6.022e23 # mol^-1
R = 8.3145 # J mol^-1 K^-1

def standardize_seq(seq, re_match=re.compile(fr"^[{''.join(BASES_WITH_U)}]+\Z")):
    """Standardizes the specified sequence and returns it, throws errors if sequence is invalid."""
    try:
        seq = str(seq).upper()
    except TypeError:
        raise TypeError(f"The sequence {seq} is not convertable to a string.")

    if re_match.match(seq):
        return seq
    else:
        if len(seq) == 0: return seq
        raise TypeError(f"The sequence {seq} contains invalid characters: {[i for i in seq if i not in BASES]}.")


def base_choices(base):
    """"""
    for new_base in BASES:
        if new_base != base:
            yield new_base


def base_choices_int(base):
    for ix in range(4):
        if BASES[ix] != base:
            yield ix



def seqlist_to_fasta(seqlist, filename):
    
    pathlib.Path(filename).parents[0].mkdir(parents=True, exist_ok=True)
    with open(filename, "w") as f:
        for i, seq in enumerate(seqlist):
            f.writelines(f">Seq{str(i).zfill(6)}")
            f.write("\n")
            f.writelines(seq)
            f.write("\n\n")


def fasta_to_seqlist(filename):
    d = []
    with open(filename) as f:
        for record in Bio.SeqIO.parse(f, 'fasta'):
            d.append(str(record.seq))
    return d


def txt_to_seqlist(filename):
    with open(filename) as f:
        d = [line.rstrip() for line in f]
    return d



def seqlist_to_fastq(seqlist, file_path, use_gzip=True):

    pathlib.Path(file_path).parents[0].mkdir(parents=True, exist_ok=True)

    if use_gzip:
        f_open = functools.partial(gzip.open, file_path, mode='wt')
    else:
        f_open = functools.partial(open, file_path, mode='w')
    

    with f_open() as f:
        f.writelines(['\n'.join([
            f"@Read{str(i).zfill(9)}",
            str(sequence),
            "+",
            "F"*len(sequence),
            '' # final newline
        ]) for i, sequence in enumerate(seqlist)])
