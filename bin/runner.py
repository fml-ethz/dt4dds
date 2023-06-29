import argparse
import subprocess

import dt4dds
dt4dds.default_logging()

parser = argparse.ArgumentParser(description='Process one experiment folder')
parser.add_argument('type', type=str, help='Type of experiment to analyse', choices=['single', 'paired', 'paired_compare', 'paired_seqanalysis', 'single_seqanalysis', 'phix', 'scafstats_single', 'scafstats_paired'])
parser.add_argument('folderpath', type=str, help='Path to the experiment folder')
parser.add_argument('-n', '--n_reads', type=int, help='Number of reads to use', default=10000)
parser.add_argument('-t', '--trim_length', type=int, help='Length to trim reads to', default=150)
parser.add_argument('-so', '--ref_offset', type=int, help='Offset for the reference sequences', default=0)
parser.add_argument('-ro', '--read_offset', type=int, help='Offset for the read sequences', default=0)
parser.add_argument('-rl', '--read_length', type=int, help='Maximum length of the reads', default=0)
parser.add_argument('-ad', '--adapter', type=str, help='Sequencing adapter', default='AGATCGGAAGAGCG')

args = parser.parse_args()


if args.type == 'paired':
    dt4dds.analysis.errormapper.paired(
        f"{args.folderpath}R1.fq.gz", 
        f"{args.folderpath}R2.fq.gz", 
        f"{args.folderpath}design_files.fasta",
        n_reads = args.n_reads,
        ref_offset = args.ref_offset,
        read_offset = args.read_offset,
        read_length = args.read_length,
        adapter = args.adapter,
    )

elif args.type == 'single':
    dt4dds.analysis.errormapper.single(
        f"{args.folderpath}R1.fq.gz", 
        f"{args.folderpath}design_files.fasta",
        n_reads = args.n_reads,
        ref_offset = args.ref_offset,
        read_offset = args.read_offset,
        read_length = args.read_length,
        adapter = args.adapter,
    )

elif args.type == 'paired_compare':
    dt4dds.analysis.errormapper.paired_cmp(
        f"{args.folderpath}R1.fq.gz", 
        f"{args.folderpath}R2.fq.gz", 
        f"{args.folderpath}design_files.fasta",
        n_reads = args.n_reads,
        ref_offset = args.ref_offset,
        read_offset = args.read_offset,
        read_length = args.read_length,
        adapter = args.adapter,
    )

elif args.type == 'paired_seqanalysis':
    dt4dds.analysis.errormapper.paired_seqanalysis(
        f"{args.folderpath}R1.fq.gz", 
        f"{args.folderpath}R2.fq.gz", 
        f"{args.folderpath}design_files.fasta",
        n_reads = args.n_reads,
        ref_offset = args.ref_offset,
        read_offset = args.read_offset,
        read_length = args.read_length,
        adapter = args.adapter,
    )

elif args.type == 'single_seqanalysis':
    dt4dds.analysis.errormapper.single_seqanalysis(
        f"{args.folderpath}R1.fq.gz", 
        f"{args.folderpath}design_files.fasta",
        n_reads = args.n_reads,
        ref_offset = args.ref_offset,
        read_offset = args.read_offset,
        read_length = args.read_length,
        adapter = args.adapter,
    )


elif args.type == 'phix':
    dt4dds.analysis.errormapper.phix(
        f"{args.folderpath}R1.fq.gz", 
        f"{args.folderpath}R2.fq.gz", 
        n_reads = args.n_reads,
        max_read_length = args.trim_length,
    )

elif args.type == 'scafstats_paired':
    cmd = f"~/.local/bin/bbmap/bbduk.sh -Xmx6g in1={args.folderpath}R1.fq.gz in2={args.folderpath}R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref={args.folderpath}design_files.fasta ordered interleaved nodisk scafstats={args.folderpath}scafstats.txt"
    output = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(cmd)
    print(output.stderr)

elif args.type == 'scafstats_single':
    cmd = f"~/.local/bin/bbmap/bbduk.sh -Xmx6g in={args.folderpath}R1.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref={args.folderpath}design_files.fasta int=t ordered nodisk scafstats={args.folderpath}scafstats.txt"
    output = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(cmd)
    print(output.stderr)
