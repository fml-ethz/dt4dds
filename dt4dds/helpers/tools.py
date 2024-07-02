import numpy as np
import pathlib
import gzip
import functools
import Bio.SeqIO


# random instance for all random operations, change via tools.set_seed(seed)
rng = np.random.default_rng(np.random.SFC64())

BASES = ["A", "C", "G", "T"]
NA = 6.022e23 # mol^-1
R = 8.3145 # J mol^-1 K^-1


def base_choices(base):
    """"""
    for new_base in BASES:
        if new_base != base:
            yield new_base


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
        f_open = functools.partial(gzip.open, file_path, mode='wt', compresslevel=1)
    else:
        f_open = functools.partial(open, file_path, mode='w')
    
    with f_open() as f:
        f.writelines(f"@Seq{i}\n{seq}\n+\n{'F'*len(seq)}\n" for i, seq in enumerate(seqlist))
