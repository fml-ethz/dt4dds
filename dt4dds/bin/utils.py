import dt4dds
import argparse
import sys
import os
import subprocess
import Bio.SeqIO
import seqfold

dt4dds.default_logging()

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# 
# FILE CONVERSION UTILITIES
# 

def txt2fasta(args):
    parser = argparse.ArgumentParser(description='Converting sequences in a text file to fasta format.')
    parser.add_argument('input', type=str, help='Path to the input file')
    parser.add_argument('output', type=str, help='Path to the output file')
    config = parser.parse_args(args)

    i = 0
    with open(config.input, 'r') as fi, open(config.output, 'w') as fo:
        for line in fi:
            line = line.strip()
            if line != "":
                i += 1
                fo.write(f">Seq{str(i).zfill(6)}\n{line}\n\n")
    logger.info(f"Converted {i} sequences from {config.input} to {config.output}.")


def fasta2txt(args):
    parser = argparse.ArgumentParser(description='Converting sequences in a fasta file to text format.')
    parser.add_argument('input', type=str, help='Path to the input file')
    parser.add_argument('output', type=str, help='Path to the output file')
    config = parser.parse_args(args)
    
    i = 0
    with open(config.output, 'w') as fo:
        for seq in Bio.SeqIO.parse(config.input, 'fasta'): 
            i += 1
            fo.write(f"{str(seq.seq)}\n")
    logger.info(f"Converted {i} sequences from {config.input} to {config.output}.")


def txt2fastq(args):
    parser = argparse.ArgumentParser(description='Converting sequences in a text file to fastq format.')
    parser.add_argument('input', type=str, help='Path to the input file')
    parser.add_argument('output', type=str, help='Path to the output file')
    config = parser.parse_args(args)

    i = 0
    with open(config.input, 'r') as fi, open(config.output, 'w') as fo:
        for line in fi:
            line = line.strip()
            if line != "":
                i += 1
                fo.write(f"@Seq{str(i).zfill(6)}\n{line}\n+\n{'F'*len(line)}\n")
    logger.info(f"Converted {i} sequences from {config.input} to {config.output}.")


def fastq2txt(args):
    parser = argparse.ArgumentParser(description='Converting sequences in a fastq file to text format.')
    parser.add_argument('input', type=str, help='Path to the input file')
    parser.add_argument('output', type=str, help='Path to the output file')
    config = parser.parse_args(args)
    
    i = 0
    with open(config.output, 'w') as fo:
        for seq in Bio.SeqIO.parse(config.input, 'fastq'): 
            i += 1
            fo.write(f"{str(seq.seq)}\n")
    logger.info(f"Converted {i} sequences from {config.input} to {config.output}.")


def fastq2fasta(args):
    parser = argparse.ArgumentParser(description='Converting sequences in a fastq file to fasta format.')
    parser.add_argument('input', type=str, help='Path to the input file')
    parser.add_argument('output', type=str, help='Path to the output file')
    config = parser.parse_args(args)
    
    i = 0
    with open(config.output, 'w') as fo:
        for seq in Bio.SeqIO.parse(config.input, 'fastq'): 
            i += 1
            fo.write(f">Seq{str(i).zfill(6)}\n{str(seq.seq)}\n\n")
    logger.info(f"Converted {i} sequences from {config.input} to {config.output}.")


def fasta2fastq(args):
    parser = argparse.ArgumentParser(description='Converting sequences in a fasta file to fastq format.')
    parser.add_argument('input', type=str, help='Path to the input file')
    parser.add_argument('output', type=str, help='Path to the output file')
    config = parser.parse_args(args)
    
    i = 0
    with open(config.output, 'w') as fo:
        for seq in Bio.SeqIO.parse(config.input, 'fasta'): 
            i += 1
            fo.write(f"@Seq{str(i).zfill(6)}\n{seq.seq}\n+\n{'F'*len(seq.seq)}\n")
    logger.info(f"Converted {i} sequences from {config.input} to {config.output}.")


# 
# FILE OPERATIONS UTILITIES
# 

def mergereads(args):
    parser = argparse.ArgumentParser(description='Merge paired reads using NGmerge.')
    parser.add_argument('input1', type=str, help='Path to the forwards read file')
    parser.add_argument('input2', type=str, help='Path to the reverse read file')
    parser.add_argument('output', type=str, help='Path to the output file')
    parser.add_argument('--min-overlap', type=int, help='Minimum overlap length', default=20)
    parser.add_argument('--min-overlap-dovetail', type=int, help='Minimum overlap of dovetailed alignments', default=20)
    parser.add_argument('--threads', type=int, help='Number of parallel threads', default=1)
    parser.add_argument('--bin', type=str, help='Path to NGmerge executable', default='~/.local/bin/NGmerge')
    config = parser.parse_args(args)

    cmd = [os.path.expanduser(config.bin), "-1", str(config.input1), "-2", str(config.input2), "-o", str(config.output), "-n", str(config.threads), "-m", str(config.min_overlap), "-e", str(config.min_overlap_dovetail), "-d", "-v"]
    logger.info(f"Running command: {' '.join(cmd)}")
    with open(f'{config.output}.ngmerge_log', 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)


def fasta2seqproperties(args):
    parser = argparse.ArgumentParser(description='Extracting sequence properties from a fasta file.')
    parser.add_argument('input', type=str, help='Path to the input file')
    parser.add_argument('output', type=str, help='Path to the output file')
    config = parser.parse_args(args)

    def get_longest_hp(seq):
        lchar = ''
        longest = 0
        for i in seq:
            if lchar == i:
                cnt += 1
            else:
                cnt = 1
            if cnt > longest:
                longest = cnt
            lchar = i
        return longest

    get_props = lambda seq: [
            seq.id,
            len(seq),
            seq.count('A')/len(seq),
            seq.count('C')/len(seq),
            seq.count('G')/len(seq),
            seq.count('T')/len(seq),
            (seq.count('G') + seq.count('C'))/len(seq),
            str(seq[0]),
            str(seq[-1]),
            get_longest_hp(str(seq.seq)),
            seqfold.dg(str(seq.seq)),
        ]
    
    i = 0
    with open(config.output, 'w') as fo:
        fo.write("id,length,A,C,G,T,GC,first,last,hp,dg\n")
        for seq in Bio.SeqIO.parse(config.input, 'fasta'): 
            i += 1
            fo.write(f"{','.join(map(str, get_props(seq)))}\n")
    logger.info(f"Analysed {i} sequences from {config.input} and saved properties to {config.output}.")







AVAILABLE_UTILITIES = {
    # File conversion utilities
    'txt2fasta': txt2fasta,
    'fasta2txt': fasta2txt,
    'txt2fastq': txt2fastq,
    'fastq2txt': fastq2txt,
    'fastq2fasta': fastq2fasta,
    'fasta2fastq': fasta2fastq,
    # File operations utilities
    'mergereads': mergereads,
    'fasta2seqproperties': fasta2seqproperties,
}


def main():
    parser = argparse.ArgumentParser(description='Running utility functions for file conversion and operations.')
    parser.add_argument('utility', type=str, help='Name of the utility to run', choices=AVAILABLE_UTILITIES.keys())
    args = parser.parse_args(sys.argv[1:2])
    logger.info(f"Running utility {args.utility}.")
    
    try:
        AVAILABLE_UTILITIES[args.utility](sys.argv[2:])
    except Exception as e:
        logger.error(f"Error running utility {args.utility}: {e}")
        sys.exit(1)